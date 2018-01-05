#include "macro.h"
subroutine internal_v()
  use openarray
  use config
  use variables
  implicit none
  type(array) dh_3d, tmpa,tmpdzz_3d
  integer :: ierr, k, i

  dh_3d=AYB(dhf_3d)

  vf =(AYB(dhb_3d)*vb-dti2*(advy+ &
       drhoy+AYB(cor_3d*dt_3d*axf_u)+     &
       grav * AYB(dt_3d)*(DYB(egf_3d+egb_3d)+ &
       DYB(e_atmos_3d)*2.d0)/2.d0 &
       -DZF(AYB(w) * AZB(v))))/dh_3d

  
  VF(kb)=0.d0

  call set(sub(dh_3d,1,':',':'),1.d0)
  call set(sub(dh_3d,':',1,':'),1.d0)

  c=AYB(km)
  !    a(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
  ! tmpa=mat_zeros
  ! call set(sub(tmpa,':',':',r(1,kbm2)),-dti2*( sub(c,':',':',r(2,kbm1))+ umol))
  ! a=tmpa/(dz_3d*dzz_3d*dh_3d*dh_3d)      !a(:,:,kbm1:kb)=0.0

  do i = 1, kbm2
     A(i) = -dti2 * (C(i+1) + umol) / (dz1(i) * dzz1(i) * DH_3D(i)**2)
  end do
  A(kbm1) = 0.d0
  A(kb) = 0.d0
  
  ! tmpdzz_3d=mat_zeros
  ! call set(sub(tmpdzz_3d,':',':',r(2,kbm1)),sub(dzz_3d,':',':',r(1,kbm2)))
  ! c=-dti2*(c+umol)/(dz_3d*tmpdzz_3d*dh_3d*dh_3d)  !c(:,:,1),c(:,:,kb)=0.0
  do i = 2, kbm2
     C(i) = -dti2 * (C(i) + umol) / (dz1(i) * dzz1(i-1) * DH_3D(i)**2)
  enddo
  
  C(1)=  0.d0
  C(kb)= 0.d0
  !call set( EE(1), A(1)/( A(1)-1.D0 ) )
  !call set( GG(1), (dti2*wvsurf/( DZ(1)*DH_3D(1))-VF(1))/(A(1)-1.d0)  )
  EE(1)=A(1)/( A(1)-1.D0 )
  GG(1)=(dti2*WVSURF_2D /( DZ(1)*DH_3D(1))-VF(1))/(A(1)-1.d0)

  !    for k=2:kbm2
  !        gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
  !        ee(:,:,k)=a(:,:,k).*gg(:,:,k);
  !        gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-vf(:,:,k)).*gg(:,:,k);
  !    end

  do k=2,kbm2
  !   call set(GG(k), 1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0))
  !   call set(EE(k), A(k) * GG(k))
  !   call set(GG(k), ( C(k) * GG(k-1) - VF(k) ) * GG(k)) 
      GG(k)= 1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0)
      EE(k)= A(k) * GG(k)
      GG(k)= ( C(k) * GG(k-1) - VF(k) ) * GG(k)
  enddo

  tps=AYB(cbc) * sqrt( AYB( AXF( UB(kbm1) ) )**2 + VB(kbm1)**2 )
  call set(sub(tps,1, ':',':'), 0.d0)
  call set(sub(tps,im,':',':'), 0.d0)

 ! call set(VF(kbm1) , ( C(kbm1)* GG(kbm2)-VF(kbm1) )/ &
 !      (tps*dti2 /(-DZ(kbm1)*sub(dh_3d,':',':',kbm1))-1.d0-C(kbm1)*( EE(kbm2)-1.d0 ) ) ) 
 VF(kbm1)= ( C(kbm1)* GG(kbm2)-VF(kbm1) )/ &
       (TPS_2D*dti2 /(-DZ(kbm1)*DH_3D(kbm1))-1.d0-C(kbm1)*( EE(kbm2)-1.d0 ) )  

  !    for k=kbm2:-1:1
  !        vf(:,:,k)=(ee(:,:,k).*vf(:,:,k+1)+gg(:,:,k));
  !    end
  do k=kbm2,1,-1
  !   call set(VF(k) , EE(k)*VF(k+1)+GG(k) )
     VF(k)= EE(k)*VF(k+1)+GG(k)
  enddo

  !    vf=vf.*dvm_3d;
  vf=vf*dvm_3d
  !    wvbot=-tps .* vf(:,:,kbm1);
  WVBOT_2D=-TPS_2D * VF(kbm1)
  !    [vf] = bcond3_v(vf, v);     

  call bcond3_v()

!  print*,"vf after bcond"
!  call disp(vf,'vf = ')

  call destroy(dh_3d, ierr);
  ! call destroy(tmpdzz_3d, ierr)
  ! call destroy(tmpa, ierr)

end subroutine
