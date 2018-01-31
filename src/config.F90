module config
  implicit none

  integer :: im, jm, kb, kl1, kl2
  integer :: kbm1, kbm2, imm1, jmm1, imm2, jmm2
  real(kind=8) :: dte, isplit, sw
  integer :: mode, nadv, nitera, nread
  real(kind=8) :: prtd1, prtd2, swtch
  integer :: iskp, jskp
  logical :: lramp
  real(kind=8) :: tprni, slmax
  integer :: ntp, nbct, nbcs, ispadv
  real(kind=8) :: alpha, tbias, sbias, tatm, satm
  integer :: iproblem
  character(len=60) :: problem
  integer :: npg
  real(kind=8) :: days, small, pi, rhoref, grav, kappa
  real(kind=8) :: z0b, cbcmin, cbcmax, horcon, umol, hmax
  real(kind=8) :: smoth, aam_init, ramp, vmaxl
  character(len=60) :: source,title, netcdf_file, time_start
  character(len=60) :: in_path
  logical :: fclim_flag, water_flag, wind_flag, heat_flag, bc_flag
  integer :: nsbdy
  
  namelist /setting/ im, jm, kb, dte, isplit, &
       kl1, kl2, iproblem, mode, nadv, problem
  namelist /setting/ nitera, sw, nread, prtd1, prtd2, days, npg
  namelist /setting/ swtch, iskp, jskp, lramp, tprni, slmax, ntp
  namelist /setting/ nbct, nbcs , ispadv, alpha , tbias, sbias, tatm, satm
  namelist /setting/ small, pi, rhoref, grav, kappa, z0b
  namelist /setting/ cbcmin, cbcmax, horcon, umol, hmax, vmaxl
  namelist /setting/ smoth, aam_init, ramp
  namelist /setting/ source, title, netcdf_file, time_start
  namelist /setting/ in_path, fclim_flag, water_flag, wind_flag, heat_flag
  namelist /setting/ nsbdy, bc_flag
  
contains

  subroutine LoadConfig()
    open(1, file='config.txt')
    read(1, setting)
    imm1 = im - 1
    imm2 = im - 2
    jmm1 = jm - 1
    jmm2 = jm - 2
    kbm1 = kb - 1
    kbm2 = kb - 2

  end subroutine

end module config
