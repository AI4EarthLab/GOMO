#include "common.h"

subroutine internal_u()
  use openarray
  use config
  use variables
  implicit none
  type(array) :: dh_3d, tmpdzz_3d, bondu
  integer :: ierr, k, i


  dh_3d=AXB( dhf_3d )
  uf=(AXB(dhb_3d)*ub-dti2*(advx + &
       drhox - AXB(cor_3d*dt_3d*ayf_v )+  &
       grav* axbdt_3d *(DXB(egf_3d+egb_3d)+ &
       DXB(e_atmos_3d)*2.d0)*0.5d0-  &
       DZF(AXB(w) * AZB(u)))) /dh_3d
  
  
  !    uf(:,:,kb)=0.e0;  bond=AXB(w) .* AZB(u); uf(im,:,:) = bond(im,:,:) ;   %add by hx
  !call set(UF(kb),0.d0)
  UF(kb)=0.d0
  bondu= AXB(w) * AZB(u)
  call set(sub(uf,im,':',':') , sub(bondu,im,':',':')) 

  call set(sub(dh_3d,1,':',':'),1.d0)
  call set(sub(dh_3d,':',1,':'),1.d0)

  !    c=AXB(km);
  c=AXB(km)
  
  !    a = create_field(zeros(im,jm,kb),gs,6);
  a = 0
  ee =mat_zeros
  gg =mat_zeros
  call grid_bind(a, gs, 6)

  ! call set(sub(a,':',':',[1,kbm2]),-dti2*( sub(c,':',':',r(2,kbm1))+ umol))
  ! a=a/(dz_3d*dzz_3d*dh_3d*dh_3d)
  
  do i = 1,kbm2
     A(i) = -dti2 * (C(i+1) + umol) / (dz1(i) * dzz1(i) * DH_3D(i)**2)
  enddo
  A(kbm1) = 0.d0
  A(kb) = 0.d0

  !  tmpdzz_3d = mat_zeros
  !!  call grid_bind(tmpdzz_3d, gs, 2)
  ! call set(sub(tmpdzz_3d,':',':',r(2,kbm1)),sub(dzz_3d,':',':',[1,kbm2]))
  ! c=-dti2*(c+umol)/(dz_3d*tmpdzz_3d*dh_3d*dh_3d)  !c(:,:,1),c(:,:,kb)=0.0

  do i = 2, kbm1
     C(i) = -dti2 * (C(i) + umol) / (dz1(i) * dzz1(i-1) * DH_3D(i)**2)
  end do
  
  C(1) = 0; C(kb) = 0
  
!   call set(EE(1),A(1)/(A(1)-1.d0) )
!   call set(GG(1),(dti2*wusurf/(DZ(1)*DH_3D(1))-UF(1))/(A(1)-1.d0))
    EE(1)=A(1)/(A(1)-1.d0)
    GG(1)=(dti2*WUSURF_2D/(DZ(1)*DH_3D(1))-UF(1))/(A(1)-1.d0)

  do k=2,kbm2
!     call set(GG(k), 1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0))
!     call set(EE(k), A(k) * GG(k))
!     call set(GG(k), (C(k) * GG(k-1) - UF(k))* GG(k))
      GG(k)=1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0)
      EE(k)=A(k) * GG(k)
      GG(K)=(C(k) * GG(k-1) - UF(k))* GG(k)
  enddo
 
  tps=AXB(cbc) * sqrt( UB(kbm1)**2 + AXB( AYF( VB(kbm1)))**2)

  
!  call set(UF(kbm1) , (C(kbm1)*GG(kbm2)-UF(kbm1))/ &
!       (tps*dti2 /(-DZ(kbm1)*sub(dh_3d,':',':',kbm1))-&
!       1.d0-C(kbm1)*( EE(kbm2)-1.d0))) 
   UF(kbm1)=(C(kbm1)*GG(kbm2)-UF(kbm1))/ &
       (TPS_2D*dti2 /(-DZ(kbm1)*DH_3D(kbm1))-&
         1.d0-C(kbm1)*( EE(kbm2)-1.d0))
  
  do k=kbm2,1,-1
     !call set(UF(k) , EE(k)*UF(k+1)+GG(k))
     UF(k)=EE(k)*UF(k+1)+GG(k)
  end do

  
  uf=uf*dum_3d;

  !wubot=-tps*UF(kbm1)
  WUBOT_2D=-TPS_2D*UF(kbm1)
  
  call bcond3_u()

  call destroy(dh_3d, ierr);
  call destroy(tmpdzz_3d, ierr)
  call destroy(bondu, ierr)
end subroutine
