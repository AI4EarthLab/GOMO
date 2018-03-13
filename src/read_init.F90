#include "common.h"

subroutine read_init(fnc)
  use openarray
  use variables
  use config
  use openarray
  use mpi
  implicit none
  character(len=*) :: fnc
  integer :: ii, ierr,i,j,k
  real*8, allocatable :: tmp_z1(:,:,:), &
       tmp_zz1(:,:,:), tmp_dz1(:,:,:), &
       tmp_dzz1(:,:,:), RR3(:,:,:)
  type(array) :: AA
  
  if(get_rank() == 0) &
       print*, "start reading init variables..."

  z       = load(fnc,'z');
  zz      = load(fnc,'zz');
  dz      = load(fnc,'dz');
  dzz     = load(fnc,'dzz');

  !covert the distributed array object to local array
  allocate(z1(kb), zz1(kb), dz1(kb), dzz1(kb))

  allocate(tmp_z1(1,1,kb), tmp_zz1(1,1,kb), &
       tmp_dz1(1,1,kb), tmp_dzz1(1,1,kb))
  
  tmp_z1   = z
  tmp_zz1  = zz
  tmp_dz1  = dz
  tmp_dzz1 = dzz

  tmpdzz=dzz
  call set(sub(tmpdzz,':',':',[2,kbm1]),sub(dzz,':',':',[1,kbm2]))

  !stencil must be correctly set!
  z1   = tmp_z1(1,1,:)
  zz1  = tmp_zz1(1,1,:)
  dz1  = tmp_dz1(1,1,:)
  dzz1 = tmp_dzz1(1,1,:)

  z_3d  = rep(z,im,jm,1);

  dx      = load(fnc,'dx')       ; !call disp(dx, 'dx = ')
  dy      = load(fnc,'dy')       ; !call disp(dy, 'dy = ')
  cor     = load(fnc,'cor')      ; !call disp(cor, 'cor = ')
  h       = load(fnc,'h')        ;
  fsm     = load(fnc,'fsm')      ;
  dum     = load(fnc,'dum')      ;
  dvm     = load(fnc,'dvm')      ;
  art     = load(fnc,'art')      ;
  aru     = load(fnc,'aru')      ; !call disp(aru, 'aru = ')
  arv     = load(fnc,'arv')      ; !call disp(arv, 'arv = ')
  rfe     = load(fnc,'rfe')      ; 
  rfw     = load(fnc,'rfw')      ;!call disp(rfw, 'rfw = ') 
  rfn     = load(fnc,'rfn')      ;!call disp(rfn, 'rfn = ')
  rfs     = load(fnc,'rfs')      ;!call disp(rfs, 'rfs = ')   
  east_e  = load(fnc,'east_e')   ;
  north_e = load(fnc,'north_e')  ;
  east_c  = load(fnc,'east_c')   ;
  north_c = load(fnc,'north_c')  ;
  east_u  = load(fnc,'east_u')   ;
  north_u = load(fnc,'north_u')  ;
  east_v  = load(fnc,'east_v')   ;
  north_v = load(fnc,'north_v')  ;
  tb      = load(fnc,'tb')       ;
  sb      = load(fnc,'sb')       ;
  tclim   = load(fnc,'tclim')    ;
  sclim   = load(fnc,'sclim')    ;
  rot     = load(fnc,'rot');
  ! uvel    = load(fnc,'uvel')     ;
  ! vvel    = load(fnc,'vvel')     ;
  vfluxf  = load(fnc,'vfluxf');
  wusurf  = load(fnc,'wusurf');
  wvsurf  = load(fnc,'wvsurf');
  e_atmos = load(fnc,'e_atmos');
  ub      = load(fnc,'ub');
  vb      = load(fnc,'vb');
  uab     = load(fnc,'uab');
  vab     = load(fnc,'vab');
  elb     = load(fnc,'elb');
  etb     = load(fnc,'etb');
  dt      = load(fnc,'dt')       ;
  uabw    = load(fnc,'uabw')     ;
  uabe    = load(fnc,'uabe')     ;
  vabs    = load(fnc,'vabs')     ;
  vabn    = load(fnc,'vabn')     ;
  els     = load(fnc,'els')      ;
  eln     = load(fnc,'eln')      ;
  ele     = load(fnc,'ele')      ;
  elw     = load(fnc,'elw')      ;
  ssurf   = load(fnc,'ssurf')    ;
  tsurf   = load(fnc,'tsurf')    ;
  tbe     = load(fnc,'tbe')      ;
  sbe     = load(fnc,'sbe')      ;
  sbw     = load(fnc,'sbw')      ;
  tbw     = load(fnc,'tbw')      ;
  tbn     = load(fnc,'tbn')      ;
  tbs     = load(fnc,'tbs')      ;
  sbn     = load(fnc,'sbn')      ;
  sbs     = load(fnc,'sbs')      ;

  wtsurf = load(fnc,'wtsurf');
  swrad  = load(fnc,'swrad');

  call dens(rho, tb, sb)

  if(fclim_flag) then
     call dens(rmean, sclim, tclim)
  else
     call set(sub(rmean, ':',':',[1,kbm1]), &
          sub(rho,':',':',[1,kbm1]))
  end if

end subroutine
