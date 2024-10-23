!>\file  mpas_smoke_config.F90
!! This file contains the configuration for MPAS-Smoke.
!
! Haiqin.Li@noaa.gov  
! 10/2024
!
module mpas_smoke_config

  use mpas_kind_types

  implicit none

  !-- constant paramters
  real(RKIND), parameter :: epsilc     = 1.e-12
  real(RKIND), parameter :: pi         = 3.1415926
  !-- aerosol module configurations
  integer :: chem_opt = 1   
  integer :: kemit = 1
  integer :: dust_opt = 1
  real(RKIND) :: dust_drylimit_factor  = 1.0
  real(RKIND) :: dust_moist_correction = 1.0
  real(RKIND) :: dust_alpha = 0.
  real(RKIND) :: dust_gamma = 0.
  integer :: seas_opt = 0   ! turn off by default
  integer :: plume_wind_eff = 0   ! turn off by default
  logical :: do_plumerise  = .true.
  integer :: addsmoke_flag = 1
  integer :: plumerisefire_frq=60
  integer :: n_dbg_lines = 3
  integer :: wetdep_ls_opt = 1
  integer :: drydep_opt  = 1
  integer :: pm_settling = 1
  integer :: nfire_types = 5
  integer :: ebb_dcycle  = 2 ! 1: read in ebb_smoke(i,24), 2: daily
  integer :: hwp_method = 2
  logical :: do_mpas_smoke = .true.
  logical :: dbg_opt       = .true.
  logical :: aero_ind_fdb  = .false.
  logical :: add_fire_heat_flux= .false.
  logical :: add_fire_moist_flux= .false.
  real(RKIND) :: wetdep_ls_alpha = .5 ! scavenging factor
  real(RKIND) :: plume_alpha = 0.05
  real(RKIND) :: sc_factor = 1.0

  ! --
  integer, parameter :: CHEM_OPT_GOCART= 1
  integer, parameter :: num_moist=2, num_chem=20, num_emis_seas=5, num_emis_dust=5

  ! -- set pointers to predefined atmospheric tracers
  integer, parameter :: ndvel = 1

  integer :: numgas = 0

  !-- tracers
   integer, parameter :: p_smoke=5
   integer, parameter :: p_dust_1=10
   integer, parameter :: p_dust_2=11
   integer, parameter :: p_dust_3=12
   integer, parameter :: p_dust_4=13
   integer, parameter :: p_dust_5=14, p_coarse_pm=14
   integer, parameter :: p_seas_1=15
   integer, parameter :: p_seas_2=16
   integer, parameter :: p_seas_3=17
   integer, parameter :: p_seas_4=18
   integer, parameter :: p_seas_5=19

   integer, parameter :: p_edust1=1,p_edust2=2,p_edust3=3,p_edust4=4,p_edust5=5

  ! -- fire options
!  integer, parameter :: num_plume_data = 1


end module
