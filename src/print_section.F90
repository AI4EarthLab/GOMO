#include "common.h"

subroutine print_section(iint,iext,vamax,imax,jmax,vtot,atot,taver,saver,eaver,tsalt)
  use openarray
  use config
  use variables
  implicit none
  real(kind=8)  eaver, taver, saver, tsalt,vtot, atot
  type(array)   darea
  real(kind=8) :: vamax, realdz
  integer :: imax, jmax, iint, iext
  integer :: ierr

  if (iint >= iswtch) then
     iprint=int(prtd2*24.0*3600.0/dti + 0.5)
  endif

!     print*, 'print_section'
!  if (mod(iint,iprint)==0 .or. vamax >= vmaxl) then
     realdz=sum(dz,3)
     darea=dx*dy*fsm
     vtot =sum(sum(realdz * darea * dt,2),1)
     atot =sum(sum(darea,2),1)

     taver=sum(sum(sum(dz*tb*darea*dt,3),2),1)/vtot
     saver=sum(sum(sum(dz*sb*darea*dt,3),2),1)/vtot
     tsalt=(saver+sbias)*vtot
     eaver=sum(sum(et*darea,2),1)/atot

!  write(60,(F21.7,F21.7,F12.7,F12.7,F12.7,F27.7)) vtot,atot, eaver,taver,saver,tsalt

     if(get_rank() == 0) then
        write(*,*) '******************************************************************************'
        write(*,'(A,F16.7,A,I8,A,i8,A,I8)') &
             " time=",time,    ";   iint=",iint, &
             ";   iext=",iext, ";   iprint=",iprint

        write(*,*) '******************************************************************************'
        write(*,'(A,F25.4,A,F22.7,A,F30.7/A,F25.10,        &
             A,F22.10,A,F30.7)') &
             '   vtot=',vtot,  ';   atot=',atot, &
             ';  eaver=',eaver, '  taver=',taver, &
             ';  saver=',saver, ';  tsalt=',tsalt
     endif
end subroutine print_section
