#include "common.h"

subroutine dens(rhoo, t_, s_)
  use openarray  
  use variables
  use config
  implicit none
  type(array), intent(inout) :: rhoo
  type(array), intent(in) :: t_, s_
  type(array) :: tr, sr, tr2, tr3, tr4, p, &
       cr, rhor, rhor1, rhor2, cr1
  integer :: ierr

  tr  = t_ + tbias;
  sr  = s_ + sbias;
  tr2 = tr * tr;
  tr3 = tr2 * tr;
  tr4 = tr3 * tr;

  p=grav * rhoref * (-zz * mat_ones *  h) * 1.d-5;

  rhor2=-0.157406d0+6.793952d-2 * tr-9.095290d-3 * tr2 &
       +1.001685d-4 * tr3 &
       -1.120083d-6 * tr4+6.536332d-9 * tr4 * tr;

  rhor1=rhor2+(0.824493d0-4.0899d-3 * tr+7.6438d-5 * tr2 &
       -8.2467d-7 * tr3 &
       +5.3875d-9 * tr4) * sr+(-5.72466d-3+1.0227d-4 * tr &
       -1.6546d-6 * tr2) * abs(sr)**1.5d0+4.8314d-4 * sr * sr;

  cr1=1449.1d0+.0821d0 * p+4.55d0 * tr-.045d0 &
       * tr2+1.34d0 * (sr-35.d0);
  
  cr=p / (cr1 * cr1)
  rhor=rhor1+1.d5 * cr * (1.d0-2.d0 * cr);
  rhoo = rhor / rhoref * fsm;
  
end subroutine
