! radiation_photolysis.F90 - Derived type containing data/routines to compute photolysis rates
!
! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

#include "ecrad_config.h"

module radiation_photolysis

  use parkind1, only : jprb
  
  implicit none
  public

  integer, parameter :: NMaxProcessNameLen = 20
  
  !---------------------------------------------------------------------
  ! This derived type contains all the data needed to calculate
  ! photolysis rates from spectral fluxes output from ecRad
  type photolysis_type

    ! Look-up table versus temperature of scaled photolysis
    ! cross-sections, so that when the matrix for a particular
    ! temperature is extracted, it can be matrix-vector multiplied by
    ! the vector of actinic fluxes in a vector of ecCKD g-points and
    ! will return a vector of photolysis rates for a set of photolysis
    ! processes. This has already been multiplied by the quantum yield
    ! so may be less than the actual absorption
    ! cross-section. Dimensioned (nproc,ng,ntemp).
    real(jprb), allocatable :: cross_section_lut(:,:,:)

    ! Store the process names as an allocatable array of fixed-length strings
    character(len=NMaxProcessNameLen), allocatable :: process_names(:)

    ! If you need to convert actinic flux in W m-2 to #photons s-1
    ! m-2, multiply by this number in J^-1, dimensioned (ng)
    real(jprb), allocatable :: photons_per_joule(:)
    
    ! Number of temperatures in look-up table
    integer :: ntemperature = 0

    ! Number of g-points in look-up table, which may be fewer than in
    ! the full ecCKD gas-optics model if only a subset are relevant
    ! for photolysis
    integer :: ng = 0

    ! Number of photolysis rates to be computed. If a particular gas
    ! has more than one process then it needs to be computed
    ! separately so is counted as a separate process. The names of
    ! processes are typically of the form <gas>_<product>, e.g. o3_o
    ! and o3_o1d to distinguish the reactions O3+hv->O2+O(3P) from
    ! O3+hv->O2+O(1D).
    integer :: nproc = 0

    ! Indices to the first and last g-point of the ecCKD model that is
    ! relevant for photolysis
    integer :: istartg = 0, iendg = 0

    ! First temperature in look-up table and spacing, in Kelvin
    real(jprb) :: temperature1, dtemperature

    ! Do we return rates on half-levels?
    logical :: do_half_level_rates = .false.

  contains
    procedure :: configure
    procedure :: calculate
    procedure :: save
    procedure :: get_absorption_from_ckd
    
  end type photolysis_type
  
  
contains

  !---------------------------------------------------------------------
  ! Configure the photolysis_type structure from a netCDF file and a
  ! list of processes to consider
  subroutine configure(this, config, file_name, processes, &
       &               iverbose, do_half_level_rates, use_ckd_for_o2)

#ifdef EASY_NETCDF_READ_MPI
    use easy_netcdf_read_mpi, only : netcdf_file    
#else
    use easy_netcdf,     only : netcdf_file
#endif
    use radiation_config,only : config_type
    use radiation_io,    only : nulout, radiation_abort
    use radiation_constants, only : PlanckConstant, SpeedOfLight
    use interpolation,   only : interpolate
    use yomhook,         only : lhook, dr_hook, jphook

    class(photolysis_type), intent(inout) :: this
    ! Name of file containing photolysis cross-sections
    type(config_type),      intent(in)    :: config
    ! Configuration structure, only used for the data directory
    character(len=*),       intent(in)    :: file_name
    ! Array of strings containing the photolysis processes that are
    ! required. These are of the form "gas[_product][:CKD]",
    ! e.g. o2_o1d:CKD, which means photolysis of O2 to produce the
    ! O(1D) oxygen radical, and to obtain the cross sections from the
    ! CKD gas optics model rather than the photolysis netCDF file.
    character(len=*),       intent(in)    :: processes(:)
    ! Verbosity level from 1 (least) to 5 (most verbose)
    integer, optional,      intent(in)    :: iverbose
    ! Do we return photolysis rates on half levels?
    logical, optional,      intent(in)    :: do_half_level_rates
    ! Do we get the cross-sections of molecular oxygen from the CKD
    ! gas optics model?  This is usualy needed because the
    ! Schumann-Runge O2 absorption lines are so fine. Default true.
    logical, optional,      intent(in)    :: use_ckd_for_o2
    
    type(netcdf_file) :: file

    ! Temperature for look-up table
    real(jprb), allocatable :: temperature(:) ! K

    ! Photolysis data for one process as read from file
    real(jprb), allocatable :: wavelength_nm(:)
    real(jprb), allocatable :: cross_section_cm2(:,:)
    real(jprb), allocatable :: quantum_yield(:,:)
    real(jprb), allocatable :: wavenumber_cm1(:) ! Reverse order

    ! Photolysis data for one process interpolated to ecCKD wavenumber
    ! grid
    real(jprb), allocatable :: wavenumber_int_cm1(:) ! cm-1
    real(jprb), allocatable :: cross_section_int_cm2(:)
    real(jprb), allocatable :: quantum_yield_int(:)
    real(jprb), allocatable :: solar_photon_flux(:) ! s-1 m-2
    real(jprb), allocatable :: photolysis_multiplier(:)
    
    ! Renormalized gpoint_fraction from ecCKD spectral definition
    ! file, dimensioned (nwavn,ng). The original gpoint_fraction sums
    ! to 1 along the wavenumber direction, whereas the new one sums to
    ! one along the g-point direction, thereby stating what fraction
    ! of the solar energy at a particular wavenumber is dealt with by
    ! each g-point.
    real(jprb), allocatable :: gpoint_fraction_renorm(:,:)

    ! Names of photolysis process and gas (e.g. o2_o1d and o2)
    character(len=NMaxProcessNameLen) :: process_name, gas_name
    ! Number of characters in process_name and gas_name
    integer :: ilenprocname, ilengasname
    
    ! Number of wavenumbers in ecCKD spectral description
    integer :: nwavn

    ! Number of wavelengths for a single process in photolysis file
    integer :: nwavl

    ! Loop indices for processes, temperatures and characters
    integer :: jproc, jt, jg, jchar
    
    ! Is the absorption cross-section of the current process
    ! temperature dependent?
    logical :: is_temperature_dependent
    
    ! Is the quantum yield of the current process temperature
    ! dependent?
    logical :: is_quantum_yield_t_dependent

    ! Do we take the absorption cross-section from the CKD file?  This
    ! is appropriate for molecular oxygen - it is triggered if the
    ! absorption cross-section is missing from the photolysis file.
    logical :: is_absorption_from_ckd
    
    integer :: iverbose_local
    logical :: use_ckd_for_o2_local
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_photolysis:configure',0,hook_handle)

    if (present(iverbose)) then
      iverbose_local = iverbose
    else
      iverbose_local = 3
    end if

    if (present(use_ckd_for_o2)) then
      use_ckd_for_o2_local = use_ckd_for_o2
    else
      use_ckd_for_o2_local = .true.
    end if
    
    associate(ckd_model => config%gas_optics_sw)
    
    if (file_name(1:1) == '/' .or. file_name(1:1) == '.') then
      ! Treat file_name as an absolute path
      call file%open(trim(file_name), iverbose=iverbose_local)
    else
      ! Assume the file is in the ecRad data directory
      call file%open(trim(config%directory_name) // '/' // trim(file_name), &
           &         iverbose=iverbose_local)
    end if

    ! Load temperature and store number, start and difference,
    ! assuming the values to be evenly spaced
    call file%get('temperature', temperature)
    this%ntemperature = size(temperature)
    this%temperature1 = temperature(1)
    this%dtemperature = (temperature(this%ntemperature)-temperature(1)) &
         &            / real(this%ntemperature,jprb)

    this%nproc = size(processes)
    allocate(this%process_names(this%nproc))
    
    ! Initially assume all g points relevant for photolysis
    this%istartg = 1
    this%iendg   = ckd_model%ng
    this%ng      = ckd_model%ng

    if (iverbose_local >= 2) then
      write(nulout,'(a,i0,a,i0,a,i0,a)') &
           &  'Setting up photolysis calculation for ', &
           &  this%nproc, ' processes and ', this%ng, &
           &  ' spectral g-points as a look-up table with ',  &
           &  this%ntemperature, ' temperatures'
    end if

    allocate(this%cross_section_lut(this%nproc,this%ng,this%ntemperature))

    ! Allocate variables on wavenumber grid
    nwavn = ckd_model%spectral_def%nwav
    allocate(wavenumber_int_cm1(nwavn))
    allocate(quantum_yield_int(nwavn))
    allocate(cross_section_int_cm2(nwavn))
    allocate(solar_photon_flux(nwavn))
    allocate(photolysis_multiplier(nwavn))
    wavenumber_int_cm1 = 0.5_jprb * (ckd_model%spectral_def%wavenumber1 &
         &                          +ckd_model%spectral_def%wavenumber2)
    wavenumber_int_cm1(1)     = ckd_model%spectral_def%wavenumber1(1)
    wavenumber_int_cm1(nwavn) = ckd_model%spectral_def%wavenumber1(nwavn)
    ! 100 converts wavenumber from units of cm-1 to m-1
    solar_photon_flux = ckd_model%spectral_def%solar_spectral_irradiance &
         &   / (PlanckConstant * SpeedOfLight * wavenumber_int_cm1 * 100.0)

    ! Renormalize gpoint_fraction
    allocate(gpoint_fraction_renorm(nwavn,this%ng))
    gpoint_fraction_renorm = ckd_model%spectral_def%gpoint_fraction &
         &  * spread(ckd_model%spectral_def%solar_spectral_irradiance, 2, this%ng)
    gpoint_fraction_renorm = gpoint_fraction_renorm &
         &  * spread(ckd_model%spectral_def%solar_irradiance &
         &           /sum(gpoint_fraction_renorm, 1), 1, nwavn)
    gpoint_fraction_renorm = gpoint_fraction_renorm &
         &  / spread(sum(gpoint_fraction_renorm, 2), 2, this%ng)

    ! Calculate conversion from actinic flux to photons s-1 m-2, if required
    allocate(this%photons_per_joule(this%ng))
    this%photons_per_joule = matmul(solar_photon_flux, gpoint_fraction_renorm) &
         &  / ckd_model%spectral_def%solar_irradiance
    if (iverbose_local >= 3) then
      write(nulout,*) '  Photons per joule in each g-point:'
      do jg = 1,this%ng
        write(nulout,'(a,i0,a,e14.8)') '    ', jg, ' ', this%photons_per_joule(jg)
      end do
    end if
    
    ! Loop over requested processes
    do jproc = 1,this%nproc
      
      is_absorption_from_ckd = .false.
      ! Find the process name, checking for a ":CKD" suffix
      ilenprocname = len_trim(processes(jproc))
      if (ilenprocname > 4) then
        if (processes(jproc)(ilenprocname-3:ilenprocname) == ':CKD') then
          process_name = processes(jproc)(1:ilenprocname-4)
          is_absorption_from_ckd = .true.
        else
          process_name = trim(processes(jproc))
        end if
      else
        process_name = trim(processes(jproc))
      end if
      this%process_names(jproc) = trim(process_name)
      
      ! Find the gas name, checking for an underscore
      ilengasname = ilenprocname
      do jchar = 2,ilenprocname
        if (process_name(jchar:jchar) == '_') then
          ilengasname = jchar - 1
          exit
        end if
      end do
      gas_name = process_name(1:ilengasname)

      if (use_ckd_for_o2_local .and. trim(gas_name) == 'o2') then
        is_absorption_from_ckd = .true.
      end if
      
      ! Load photolysis data for this process
      call file%get(trim(process_name) // "_wavelength", wavelength_nm)
      call file%get(trim(process_name) // "_quantum_yield", quantum_yield)

      if (.not. is_absorption_from_ckd) then
        call file%get(trim(process_name) // "_cross_section", cross_section_cm2)
        if (size(cross_section_cm2,2) == this%ntemperature) then
          is_temperature_dependent = .true.
        else
          is_temperature_dependent = .false.
        end if
      else
        is_temperature_dependent = .true.
      end if

      ! Interpolate on to ecCKD wavenumber grid
      nwavl = size(wavelength_nm)

      allocate(wavenumber_cm1(nwavl))
      wavenumber_cm1 = 1.0e7_jprb/wavelength_nm(nwavl:1:-1)

      if (size(quantum_yield,2) == this%ntemperature) then
        is_quantum_yield_t_dependent = .true.
      else
        is_quantum_yield_t_dependent = .false.
      end if
      
      if (iverbose_local >= 2) then
        if (is_temperature_dependent .or. is_quantum_yield_t_dependent) then
          write(nulout,'(a,a,a,f0.1,a,f0.1,a)') '  Temperature-dependent photolysis "', &
               &  trim(process_name), &
               &  '" sensitive to wavenumbers ', wavenumber_cm1(1), '-', &
               &  wavenumber_cm1(nwavl), ' cm-1'
        else
          write(nulout,'(a,a,a,f0.1,a,f0.1,a)') '  Temperature-independent photolysis "', &
               &  trim(process_name), &
               &  '" sensitive to wavenumbers ', wavenumber_cm1(1), '-', &
               &  wavenumber_cm1(nwavl), ' cm-1'
        end if
      end if
      
      if (wavenumber_int_cm1(1) > wavenumber_cm1(1) &
           &  .or. wavenumber_int_cm1(nwavn) < wavenumber_cm1(nwavl)) then
        if (iverbose >= 1) then
          write(nulout,'(a,a,a,f0.1,a,f0.1,a)') '    Warning: photolysis of "', &
               &  trim(process_name), &
               &  '" sensitive to wavenumbers out of available range ', &
               &  wavenumber_int_cm1(1), '-', wavenumber_int_cm1(nwavn), ' cm-1'
        end if
      end if

      ! Need to cope with the case of either of the absorption
      ! cross-section or the quantum yield being temperature dependent
      if (is_absorption_from_ckd) then
        call this%get_absorption_from_ckd(trim(gas_name), ckd_model, &
             &                            this%cross_section_lut(jproc,:,:))
        do jt = 1,this%ntemperature
          if (is_quantum_yield_t_dependent) then
            call interpolate(wavenumber_cm1, quantum_yield(nwavl:1:-1,jt), &
                 &           wavenumber_int_cm1, quantum_yield_int, 0.0_jprb)
          else
            call interpolate(wavenumber_cm1, quantum_yield(nwavl:1:-1,1), &
                 &           wavenumber_int_cm1, quantum_yield_int, 0.0_jprb)
          end if
          photolysis_multiplier = quantum_yield_int * solar_photon_flux;
          this%cross_section_lut(jproc,:,jt) = this%cross_section_lut(jproc,:,jt) &
               &  * matmul(photolysis_multiplier, gpoint_fraction_renorm) &
               &  / ckd_model%spectral_def%solar_irradiance
        end do
      else if (is_temperature_dependent .or. is_quantum_yield_t_dependent) then
        do jt = 1,this%ntemperature
          ! Interpolate to wavenumber grid
          if (is_temperature_dependent) then
            call interpolate(wavenumber_cm1, cross_section_cm2(nwavl:1:-1,jt), &
                 &           wavenumber_int_cm1, cross_section_int_cm2, 0.0_jprb)
          else
            call interpolate(wavenumber_cm1, cross_section_cm2(nwavl:1:-1,1), &
                 &           wavenumber_int_cm1, cross_section_int_cm2, 0.0_jprb)
          end if
          if (is_quantum_yield_t_dependent) then
            call interpolate(wavenumber_cm1, quantum_yield(nwavl:1:-1,jt), &
                 &           wavenumber_int_cm1, quantum_yield_int, 0.0_jprb)
          else
            call interpolate(wavenumber_cm1, quantum_yield(nwavl:1:-1,1), &
                 &           wavenumber_int_cm1, quantum_yield_int, 0.0_jprb)
          end if
          ! 1e-4 converts cross section from cm2 to m2
          photolysis_multiplier = 1.0e-4_jprb * cross_section_int_cm2 * quantum_yield_int &
               &                * solar_photon_flux;
          this%cross_section_lut(jproc,:,jt) &
               &  = matmul(photolysis_multiplier, gpoint_fraction_renorm) &
               &  / ckd_model%spectral_def%solar_irradiance
        end do
      else
        ! No temperature dependence in either absorption cross section
        ! or quantum yield
        call interpolate(wavenumber_cm1, cross_section_cm2(nwavl:1:-1,1), &
             &           wavenumber_int_cm1, cross_section_int_cm2, 0.0_jprb)
        call interpolate(wavenumber_cm1, quantum_yield(nwavl:1:-1,1), &
             &           wavenumber_int_cm1, quantum_yield_int, 0.0_jprb)
        photolysis_multiplier = 1.0e-4_jprb * cross_section_int_cm2 * quantum_yield_int &
             &                * solar_photon_flux;
        this%cross_section_lut(jproc,:,1) &
             &  = matmul(photolysis_multiplier, gpoint_fraction_renorm) &
             &  / ckd_model%spectral_def%solar_irradiance
        this%cross_section_lut(jproc,:,2:this%ntemperature) &
             &  = spread(this%cross_section_lut(jproc,:,1),2,this%ntemperature-1)
      end if
      
      deallocate(wavelength_nm)
      deallocate(wavenumber_cm1)
      deallocate(quantum_yield)
      if (allocated(cross_section_cm2)) then
        deallocate(cross_section_cm2)
      end if
    end do
    
    if (present(do_half_level_rates)) then
      this%do_half_level_rates = do_half_level_rates
    else
      this%do_half_level_rates = .false.
    end if
 
    call file%close()

    end associate
    
    if (lhook) call dr_hook('radiation_photolysis:configure',1,hook_handle)

  end subroutine configure

  
  !---------------------------------------------------------------------
  ! Calculate photolysis rates from spectral fluxes
  subroutine calculate(this, icol, mu0, temperature_hl, flux, rates, ilay1, ilay2)

    use radiation_flux,  only : flux_type
    use radiation_io,    only : nulerr, radiation_abort
    use yomhook,         only : lhook, dr_hook, jphook

    class(photolysis_type), intent(in)  :: this
    ! Column number
    integer,                     intent(in)  :: icol
    ! Cosine of solar zenith angle
    real(jprb),                  intent(in)  :: mu0 
    ! Half-level temperature (K)
    real(jprb),                  intent(in)  :: temperature_hl(:) ! (nlay+1)
    ! Structure containing spectral fluxes from ecRad
    type(flux_type),             intent(in)  :: flux
    ! Output photodissociation rates (s-1)
    real(jprb),                  intent(out) :: rates(:,:) ! (nproc,nlay)
    ! Optional range of layers to process    
    integer, optional,           intent(in)  :: ilay1, ilay2

    ! Actinic flux at a half-level
    real(jprb) :: actinic_flux(this%ng)

    ! Local cross-sections
    real(jprb) :: cross_section(this%nproc,this%ng)

    ! Rates at half-levels
    real(jprb), allocatable :: rates_hl(:,:)
    
    integer :: il1, il2

    integer :: nlay

    integer :: jlev, jlay, jproc

    ! Index
    integer :: itemp
    real(jprb) :: wtemp
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_photolysis:calculate',0,hook_handle)

    ! Checks
    if (.not. allocated(flux%sw_dn_direct_band)) then
      write(nulerr, '(a)') '*** Error: spectral shortwave radiative fluxes not output by ecRad for photolysis'
      call radiation_abort('Radiation configuration error')
    else if (size(flux%sw_dn_direct_band,1) /= this%ng) then
      write(nulerr,'(a,i0,a,i0)') '*** Error: photolysis expects spectral fluxes at ', this%ng, &
           &  ' g-points, got ', size(flux%sw_dn_direct_band,1)
      call radiation_abort('Radiation configuration error')
    end if
      
    ! Check the sun is above the horizon
    if (mu0 > 0.0_jprb) then
    
      if (present(ilay1)) then
        il1 = ilay1
      else
        il1 = 1
      end if

      if (present(ilay2)) then
        il2 = ilay2
      else
        il2 = size(flux%sw_dn_band,3)-1
      end if

      nlay = il2-il1+1

      allocate(rates_hl(this%nproc,nlay+1))
      
      ! Loop over half-levels
      do jlev = il1,il2+1
        ! Assume the diffuse fluxes are isotropic in each hemisphere so
        ! the diffuse fluxes into a horizontal plane are multiplied by 2
        ! and the direct flux into a horizontal plane is scaled by 1/mu0
        actinic_flux = flux%sw_dn_direct_band(this%istartg:this%iendg,icol,jlev) &
             &            * (1.0_jprb/mu0 - 2.0_jprb) &
             &  + 2.0_jprb * (flux%sw_dn_band(this%istartg:this%iendg,icol,jlev) &
             &               +flux%sw_up_band(this%istartg:this%iendg,icol,jlev))
        ! Interpolation points and weights
        wtemp = 1.0_jprb + (temperature_hl(jlev) - this%temperature1) / this%dtemperature
        if (wtemp < 1.0_jprb) then
          itemp = 1
          wtemp = 0.0_jprb
        else if (wtemp >= this%ntemperature) then
          itemp = this%ntemperature-1
          wtemp = 1.0_jprb
        else
          itemp = int(wtemp)
          wtemp = wtemp - itemp
        end if
        ! Interpolate cross sections
        cross_section = (1.0_jprb - wtemp) * this%cross_section_lut(:,:,itemp) &
             &        +             wtemp  * this%cross_section_lut(:,:,itemp+1)
        
        rates_hl(:,jlev) = matmul(cross_section, actinic_flux)
      end do

      if (this%do_half_level_rates) then
        ! Copy the half-level rates to the output array
        rates = rates_hl
      else
        ! Loop over full levels
        do jlay = il1,il2
          ! Assume the photolysis rates for each process vary exponentially within each layer
          do jproc = 1,this%nproc
            if (rates_hl(jproc,jlay+1) > 0.99_jprb * rates_hl(jproc,jlay)) then
              ! Optically thin layer: very small vertical variation of
              ! rates so take average of values at top and bottom
              ! of layer
              rates(jproc,jlay) = 0.5_jprb * (rates_hl(jproc,jlay) + rates_hl(jproc,jlay+1))
            else if (rates_hl(jproc,jlay) <= 0.0_jprb) then
              ! No flux
              rates(jproc,jlay) = 0.0_jprb
            else
              ! Assume an exponential variation of actinic flux through
              ! the layer and calculate the layer-mean value
              rates(jproc,jlay) = (rates_hl(jproc,jlay+1) - rates_hl(jproc,jlay)) &
                   &      / log(max(rates_hl(jproc,jlay+1)/rates_hl(jproc,jlay),tiny(1.0_jprb)))
            end if
          end do
        end do
      end if
    else
      ! Sun below horizon
      rates(:,:) = 0.0_jprb
    end if
    
    if (lhook) call dr_hook('radiation_photolysis:calculate',1,hook_handle)
    
  end subroutine calculate

  
  !---------------------------------------------------------------------
  ! Save computed photolysis rates to a netCDF file
  subroutine save(this, file_name, rates, iverbose)

    use easy_netcdf,     only : netcdf_file
    !use radiation_io,    only : nulout, nulerr, radiation_abort
    use yomhook,         only : lhook, dr_hook, jphook

    class(photolysis_type), intent(inout) :: this
    ! Name of file containing photolysis cross-sections
    character(len=*),            intent(in)    :: file_name
    ! Photolysis rates (s-1) dimensioned (nproc,nlay,ncol)
    real(kind=jprb), allocatable :: rates(:,:,:)
    ! Verbosity level from 1 (least) to 5 (most verbose)
    integer, optional,     intent(in)    :: iverbose

    ! Object for output NetCDF file
    type(netcdf_file) :: out_file

    integer :: nlev, ncol
    integer :: jproc
    integer :: i_local_verbose

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_photolysis:save',0,hook_handle)

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = 3
    end if

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose)

    ! Define dimensions
    ncol = size(rates,3)
    nlev = size(rates,2)
    
    call out_file%define_dimension("column", ncol)
    if (this%do_half_level_rates) then
      call out_file%define_dimension("half_level",  nlev)
    else
      call out_file%define_dimension("level",  nlev)
    end if

    ! Put global attributes
    call out_file%put_global_attributes( &
         &   title_str="Photolysis rates computed from the ecRad offline radiation model", &
         &   references_str="Hogan, R. J., and A. Bozzo, 2018: A flexible and efficient radiation " &
         &   //"scheme for the ECMWF model. J. Adv. Modeling Earth Sys., 10, 1990–2008", &
         &   source_str="ecRad offline radiation model")

    do jproc = 1,this%nproc
      if (this%do_half_level_rates) then
        call out_file%define_variable(trim(this%process_names(jproc)) // "_photolysis_rate_hl", &
             &  long_name=trim(this%process_names(jproc)) // " photolysis rate at half-levels", units_str="s-1", &
             &  dim2_name="column", dim1_name="half_level")
      else
        call out_file%define_variable(trim(this%process_names(jproc)) // "_photolysis_rate", &
             &  long_name=trim(this%process_names(jproc)) // " photolysis rate", units_str="s-1", &
             &  dim2_name="column", dim1_name="level")
      end if
    end do

    do jproc = 1,this%nproc
      if (this%do_half_level_rates) then
        call out_file%put(trim(this%process_names(jproc)) // "_photolysis_rate_hl", rates(jproc,:,:))
      else
        call out_file%put(trim(this%process_names(jproc)) // "_photolysis_rate", rates(jproc,:,:))
      end if
    end do

    call out_file%close()

    if (lhook) call dr_hook('radiation_photolysis:save',1,hook_handle)

  end subroutine save

  !---------------------------------------------------------------------
  ! Get the absorption coefficients for each g-point from the CKD
  ! model for a particular gas
  subroutine get_absorption_from_ckd(this, gas_name, ckd_model, cross_section)

    use radiation_io,    only : nulout, nulerr, radiation_abort
    use yomhook,         only : lhook, dr_hook, jphook
    use radiation_gas_constants, only : GasLowerCaseName, NMaxGases
    use radiation_constants, only : AvogadroConstant
    use radiation_ecckd, only : ckd_model_type
    use radiation_ecckd_gas

    ! Photolysis calculations do not account for the pressure of the
    ! gas, so we have to choose which pressure to read from the CKD
    ! look-up table.  Since this routine is likely to only be called
    ! for molecular oxygen which is important only near the model top,
    ! we take the second pressure in the CKD look-up table from the
    ! top.
    !integer, parameter :: IRefPressure = 2
    integer, parameter :: IRefPressure = 10
    
    class(photolysis_type), intent(inout)      :: this
    character(len=*),       intent(in)         :: gas_name
    type(ckd_model_type),   intent(in), target :: ckd_model
    real(jprb),             intent(out)        :: cross_section(:,:) ! (ng,ntemp)

    ! Molar absorption, dimensioned (ng,ntemp)
    real(jprb), pointer :: molar_abs(:,:)

    ! Index to and weight of first interpolation point; second is i1+1, (1.0-weight1)
    integer :: i1
    real(jprb) :: weight1, temperature
    
    integer :: i_gas_code, jgas, jtemp
    integer :: i_target_gas, i_composite_gas    
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_photolysis:get_absorption_from_ckd',0,hook_handle)

    ! First convert gas_name to a numerical code
    i_gas_code = 0
    do jgas = 1, NMaxGases
      if (gas_name == trim(GasLowerCaseName(jgas))) then
        i_gas_code = jgas
        exit
      end if
    end do
    
    if (i_gas_code == 0) then
      write(nulerr,'(a,a,a)') '*** Error: absorption coefficients of gas "', &
           &  gas_name, '" are not available in ecCKD to use for photolysis'
      call radiation_abort('Radiation configuration error')
    end if
    
    ! Now find it in the list of CKD gases
    i_target_gas    = 0
    i_composite_gas = 0
    do jgas = 1, ckd_model%ngas
      associate (single_gas => ckd_model%single_gas(jgas))
        if (single_gas%i_gas_code == i_gas_code) then
          i_target_gas = jgas
        else if (single_gas%i_gas_code == 0) then
          i_composite_gas = jgas
        end if
      end associate
    end do

    if (i_target_gas == 0) then    
      if (i_composite_gas == 0) then
        write(nulerr,'(a,a,a)') '*** Error: could not find absorption coefficients of gas "', &
             &  gas_name, '" for photolysis in current ecCKD gas-optics scheme'
        call radiation_abort('Radiation configuration error')
      elseif (i_gas_code /= IO2) then
        write(nulerr,'(a,a,a)') '*** Error: could not find absorption coefficients of gas "', &
             &  gas_name, '" for photolysis in current ecCKD gas-optics scheme, and only O2 can be taken from the "composite" gas'
        call radiation_abort('Radiation configuration error')
      else
        write(nulout,'(a)') 'Warning: assuming ecCKD "composite" gas includes molecular oxygen for photolysis cross sections'
        i_target_gas = i_composite_gas
      end if
    else
      write(nulout,'(a,a,a)') '  Reading ', trim(gas_name), ' photolysis cross sections from ecCKD'
    end if

    ! Extract absorption cross sections
    nullify(molar_abs)
    associate (single_gas => ckd_model%single_gas(i_target_gas))
      if (single_gas%i_conc_dependence == IConcDependenceRelativeLinear) then
        write(nulerr,'(a,a,a)') '*** Error: cannot extract photolysis cross section of "', gas_name, &
             &  '" from CKD gas with a relative-linear concentration dependence'
        call radiation_abort('Radiation configuration error')
      else if  (single_gas%i_conc_dependence == IConcDependenceLUT) then
        ! Interpolate from the molar absorption at the lowest
        ! concentration level and the second highest pressure level
        molar_abs => single_gas%molar_abs_conc(:,IRefPressure,:,1)
      else
        ! Interpolate from the molar absorption at the second
        ! highest pressure level
        molar_abs => single_gas%molar_abs(:,IRefPressure,:)
      end if

      do jtemp = 1,this%ntemperature
        temperature = this%temperature1 + this%dtemperature*(jtemp-1)
        i1 = 1 + floor((temperature - ckd_model%temperature1(IRefPressure)) / ckd_model%d_temperature)
        if (i1 <= 0) then
          i1 = 1
          weight1 = 1.0_jprb
        else if (i1 >= ckd_model%ntemp) then
          i1 = ckd_model%ntemp-1
          weight1 = 0.0_jprb
        else
          weight1 = (ckd_model%temperature1(IRefPressure) + i1*ckd_model%d_temperature - temperature) / ckd_model%d_temperature
        end if
        cross_section(:,jtemp) = weight1*molar_abs(:,i1) + (1.0_jprb-weight1)*molar_abs(:,i1+1)
      end do

      nullify(molar_abs)
          
      ! Convert from molar cross section to molecular cross section
      if (i_target_gas == i_composite_gas) then
        ! Assume we are dealing with oxygen as part of composite gas 
        cross_section = cross_section / (O2NominalMoleFraction * AvogadroConstant)
      else
        cross_section = cross_section / AvogadroConstant
      end if
    end associate

    if (lhook) call dr_hook('radiation_photolysis:get_absorption_from_ckd',1,hook_handle)

  end subroutine get_absorption_from_ckd
  
end module radiation_photolysis

