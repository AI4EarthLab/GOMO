#include "common.h"

! function [drhox,drhoy] = baropg(rho, rmean, dt, ramp)
! Calculates  baroclinic pressure gradient.
subroutine baropg()
  use openarray
  use variables
  use config
  implicit none
  integer:: ierr

  zz = make_psudo3d(zz)
  dz = make_psudo3d(dz)
  dum = make_psudo3d(dum)
  
  ! drhox=ramp * grav * AXB(dt)  * csum(-DZB(zz) &
  !      * DXB(AZB(rho - rmean)) * AXB(dt) &
  !      + DXB(dt) * DZB(AXB(rho - rmean)) &
  !      * AZB(zz) * AZB(dz) , 3) * dum;

  print*, "ramp=", ramp
  ! drhox=ramp*AXB(dt)*mat_ones;

  drhox=AXB(dt)*mat_ones*ramp;
  
  call disp(drhox, 'drhox = ')
  
  drhoy=ramp * grav * AYB(dt)  * csum(-DZB(zz) &
       * DYB(AZB(rho - rmean)) * AYB(dt)  &
       + DYB(dt) * DZB(AYB(rho - rmean))  &
       * AZB(zz) * AZB(dz) , 3) * dvm;
  
end subroutine

subroutine baropg_mcc()
  use openarray
  use variables
  use config

  implicit none
  integer:: ierr
  type(array) :: d4, ddx, drho, rhou, tmp_drho
  type(array) :: tmpd4, tmpddx, tmpdrho, tmprhou
  
  ! % calculate  baroclinic pressure gradient
  ! 4 th order correction terms, following McCalpin
  ! global im kb dum dvm dx dy dz grav zz gs;

  rho=rho-rmean;

  tmpdrho=  DXB(rho) * dum; tmprhou=  AXB(rho) *dum;
  tmpddx=   DXB(d)* dum; tmpd4= AXB(d)*dum;

  !tmpdrho=create_field(tmpdrho),gs,3);     
  call grid_bind(tmp_drho, 3)
  
  drho=tmpdrho -(DXB(DXF(DXB(rho)*dum)) / AXB(dx)&
       -2.0*DXB(rho) *(1-dum))/24.0;
     
  rhou=tmprhou-DXB(2.0 * AXF(DXB(rho)*dum))/16.0;
  
  ddx =tmpddx -(DXB(DXF(DXB(d) * dum)) &
       / AXB(dx)-2.0*DXB(d) *(1-dum)) /24.0;
   
  d4  =tmpd4-DXB(2.0*AXF(DXB(d)*dum))/16.0;
  
  !drho(1:2,:,:)   =tmpdrho(1:2,:,:);       
  call set(sub(drho, [1,2], ':',':'), sub(tmpdrho, [1,2],':',':'))
  !rhou(1:2,:,:)  =tmprhou(1:2,:,:);  
  call set(sub(rhou, [1,2], ':',':'), sub(tmprhou, [1,2],':',':'))
  
  !ddx(1:2,:,:)    =tmpddx(1:2,:,:);
  call set(sub(ddx, [1,2], ':',':'), sub(tmpddx, [1,2],':',':'))
  
  !d4(1:2,:,:)    =tmpd4(1:2,:,:);
  call set(sub(d4, [1,2], ':',':'), sub(tmpd4, [1,2],':',':'))
  
  drhox=AXB(dt) * ramp * grav &
       * csum(-DZB(zz)* d4 * AZB(drho) &
       +  AZB(zz) * ddx * DZB(rhou)*AZB(dz), 3)* dum;
  
  !drhox(im,:,:)=0.0;
  call set(sub(drhox, im, ':',':'), 0)
  
  tmpdrho=   DYB(rho) * dvm;
  tmprhou=   AYB(rho) * dvm;
  
  tmpddx=    DYB(d)*dvm;
  tmpd4=     AYB(d)*dvm;
  
  !tmpdrho=create_field(tmpdrho),gs,3);
  call grid_bind(tmpdrho, 3)

  
  drho=tmpdrho-(DYB(DYF(DYB(rho) * dvm)) &
       / AYB(dy) - 2.0*DYB(rho)*(1-dvm))/24.0;
  
  rhou=tmprhou-DYB(2.0 * AYF(DYB(rho)  * dvm))/16.0;
  
  ddx =tmpddx-(DYB(DYF(DYB(d) * dvm)) &
       / AYB(dy)-2.0*DYB(d)*(1-dvm))/24.0;
  
  !d4  =tmpd4-DYB(2.0*AYF(DYB(d))*dvm)))/16.0;
  d4  =tmpd4  -DYB(2.0 * AYF(DYB(d ) * dvm))/16.0
  
  !drho(:,1:2,:)   =tmpdrho(:,1:2,:);
  call set(sub(drho, ':',[1,2],':'), sub(tmpdrho, ':',[1,2],':'))
  
  !rhou(:,1:2,:)  =tmprhou(:,1:2,:);
  call set(sub(rhou, ':',[1,2],':'), sub(tmprhou, ':',[1,2],':'))
  
  !ddx(:,1:2,:)    =tmpddx(:,1:2,:);
  call set(sub(ddx, ':',[1,2],':'), sub(tmpddx, ':',[1,2],':'))
  
  !d4(:,1:2,:)    =tmpd4(:,1:2,:);        
  call set(sub(d4, ':',[1,2],':'), sub(d4, ':',[1,2],':'))
  
  drhoy=AYB(dt) * ramp * grav &
       * csum(-DZB(zz)* d4 * AZB(drho)  &
       + AZB(zz) * ddx * DZB(rhou)*AZB(dz), 3)* dvm;
  !drhoy(im,:,:)=0.0;
  call set(sub(drhoy, im, ':', ':'), 0)

end subroutine
