#include "common.h"

subroutine print_section(iint,iext,vamax,imax,jmax,vtot,atot,taver,saver,eaver,tsalt)
  use openarray
  use config
  use variables
  implicit none
  real(kind=8)  eaver, taver, saver, tsalt
  type(array) ::  tmp_dvol, darea
  real(kind=8) :: vtot, atot
  real(kind=8) :: vamax
  integer :: imax, jmax, iint, iext
  integer :: ierr
  
  if (iint >= iswtch) then
     iprint=int(prtd2*24.0*3600.0/dti + 0.5)
  endif

!  if (mod(iint,iprint)==0 .or. vamax >= vmaxl) then
     tmp_dvol=dx_3d * dy_3d * fsm_3d * dt_3d * dz_3d
     vtot=sum(tmp_dvol)
     darea=dx*dy*fsm
     atot=sum(darea)

     taver=sum(tb*tmp_dvol)/vtot
     saver=sum(sb*tmp_dvol)/vtot
     tsalt=(saver+sbias)*vtot
     eaver=sum(et*darea)/atot

!  write(60,(F21.7,F21.7,F12.7,F12.7,F12.7,F27.7)) vtot,atot, eaver,taver,saver,tsalt

  if (mod(iint,iprint)==0 .or. vamax >= vmaxl) then
     if(get_rank() == 0) then 
        write(*,*) '******************************************************************************'
        write(*,'(A,F16.7,A,I8,A,i8,A,I8)') &
             " time=",time,    ";   iint=",iint, &
             ";   iext=",iext, ";   iprint=",iprint

        write(*,*) '******************************************************************************'
        write(*,'(A,F25.4,A,F22.7,A,F27.7/A,F25.10,		&
             A,F22.10,A,F27.7)') &
             '   vtot=',vtot,  ';   atot=',atot, &
             ';  eaver=',eaver, '  taver=',taver, &
             ';  saver=',saver, ';  tsalt=',tsalt
!        call disp(va,'va= ')
     endif

     if (vamax > vmaxl) then
        if(get_rank() == 0) then
           write(*,*) '*********************************************************'
           write(*,'(A,F16.7,A,I8,A,I8,A,I8)')"time=",time,";   iint=",iint, &
               ";   iext=",iext,";   iprint=",iprint
           write(*,*) '************************************************'
           write(*,*) '************ abnormal job end ******************' 
           write(*,*) '************* user terminated ******************'
           write(*,*) '************************************************'
           write(*,'(A,E12.3, A,I5,A,I5)') &
                'vamax=', vamax, ';    imax=', imax, ';    jmax=', jmax
!           call disp(va,'va= ')
!           call disp(vaf,'vaf= ')
        endif
        stop
     endif
  endif

  call destroy(tmp_dvol, ierr)
  call destroy(darea, ierr)
end subroutine print_section
