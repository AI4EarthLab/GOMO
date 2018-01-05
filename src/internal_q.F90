#include "macro.h"
!#include "common.h"

subroutine internal_q()
  use openarray
  use config
  use variables
  implicit none
  real*8 a1,b1,a2,b2,c1,e1,e2,sef,cbcnst,surfl,shiw
  type(array) :: utau2,l0,p,cc,boygr,gh, stf
  type(array) :: filter
  type(array) :: dh, kn, sh, sm, kappa_l0
  integer :: ierr, k, i
  
  ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
  dh = h + etf;

  q2f= (q2b*dhb-dti2*(-DZB(AZF(w*q2)) &
       + DXF(AXB(q2)* AZB(u) *AXB(dt) &
       -AZB(AXB(aam))*AXB(h)* dum *DXB(q2b)) &
       +DYF(AYB(q2)*  AZB(v)* AYB(dt) &
       -AZB(AYB(aam))*AYB(h)* dvm *DYB(q2b))))/dhf

  q2lf= (q2lb*dhb-dti2*(-DZB(AZF(w*q2l)) &
       +DXF(AXB(q2l)* AZB(u) *AXB(dt) &
       -AZB(AXB(aam))*AXB(h)* dum *DXB(q2lb)) &
       +DYF(AYB(q2l)* AZB(v)* AYB(dt) &
       -AZB(AYB(aam))*AYB(h)* dvm *DYB(q2lb))))/dhf

  call tic("q2f")
   q2f= DZB(AZF(w*q2)) &
        + DXF(AXB(q2)*AXB(dt)*AZB(u) &
        - AZB(AXB(aam))*AXB(h)*DXB(q2b)*dum) &
        + DYF(AYB(q2)*AYB(dt)*AZB(v) &
        - AZB(AYB(aam))*AYB(h)*DYB(q2b)*dvm)
  call toc("q2f")

  a1=0.92d0;      b1=16.6d0;  a2=0.74d0;
  b2=10.1d0;      c1=0.08d0
  e1=1.8d0;       e2=1.33d0;  sef=1.d0;
  cbcnst=100.d0;  surfl=2.d5; shiw=0.d0

!  a=mat_zeros;; c=mat_zeros
  ee=mat_zeros; gg=mat_zeros; l0=mat_zeros;

  !call set(sub(a,':',':',r(2,kbm1)), sub(dzz,':',':',r(1,kbm2)))
  do i = 2, kbm1
     A(i) = dzz1(i-1)
  end do

  a=-dti2*(AZF(kq)+umol)/(a*dz*dhf*dhf)

  ! AZF_kq = AZF(kq)
  ! do i = 2, kbm1
  !    A(i) = -dti2 * (AZF_kq(:,:,i-1) + umol) / &
  !         (A(i) * dz1(i) * dhf%data(:,:,0)**2)
  ! end do
  
  ! a(1)=a(kbm1:kb)=0
  call set(A(1),  0.d0)
  call set(A(kb), 0.d0)
  
  ! A(1)= 0.d0
  ! A(kb)=0.d0

  ! call set(sub(c,':',':',r(2,kbm1)), &
  !      sub(dzz,':',':',r(1,kbm2))*sub(dz,':',':',r(1,kbm2)))
  ! do i = 2, kbm1
  !    C(i) = dzz1(i-1) * dz1(i-1)
  ! end do

  call set(sub(c,':',':',[2,kbm1]), &
       sub(dzz,':',':',[1,kbm2])*sub(dz,':',':',[1,kbm2]))
  
  call set(c, shift(dzz * dz, 0,0,1));
  
  !c=-dti2*(AZB(kq)+umol)/(c*dhf*dhf)
  AZB_kq = AZB(kq)
  do i = 2, kbm1
     C(i) = -dti2 * (AZB_kq(:,:,i-1) + umol)/&
          (C(i) * dhf%data(:,:,0)**2)
  end do
  
  call set(C(1), 0.d0)
  call set(C(kb),0.d0)
  ! C(1)= 0.d0
  ! C(kb)=0.d0

  utau2 = sqrt(AXF(wusurf)**2 +AYF(wvsurf)**2)
  !    gg(:,:,1)=(15.8*cbcnst)^(2./3.)*utau2;

  call set(GG(1),(15.8d0*cbcnst)**(2.d0/3.d0)*utau2)
  !GG(1)=(15.8d0*cbcnst)**(2.d0/3.d0)*UTAU2

  !l0 = surfl/grav*utau2
  ! do k=1,kb
  !    L0(k)=surfl/grav*UTAU2
  ! enddo

  l0 = rep(surfl/grav*utau2, 1, 1, kb)
  
  ! do k=1,kb
  !    call set(L0(k), surfl/grav*utau2)
  ! enddo
  
  !  l0=surfl*utau2/grav
  !    q2f(:,:,kb) = sqrt(AXF(wubot).^2 +AYF(wvbot).^2) .* (16.6e0^(2.e0/3.e0))*sef;

  call set(sub(q2f,':',':',kb), &
       sqrt(AXF(wubot)**2 +AYF(wvbot)**2)*(16.6d0**(2.d0/3.d0))*sef)


  p=grav*rhoref*(-zz * h)*1.d-4

  cc=1449.10d0+0.00821d0*p+4.55d0*(t+tbias) &
       -0.045d0*(t+tbias)**2 +1.34d0*(s+sbias-35.0d0)
  
  cc=cc/sqrt((1.d0-0.01642d0*p/cc)*(1.d0-0.4d0*p/cc**2))

  call set(sub(cc,':',':',kb), 0.0)

  !q2b =abs(q2b);  q2lb=abs(q2lb);
  q2b =abs(q2b);;  q2lb=abs(q2lb)
  !boygr=-grav* DZB(rho)./h + DIVISION(grav^2 , AZB(cc.^2));  boygr(:,:,1)=0.e0;
  boygr=-grav*DZB(rho)/h + grav**2/AZB(cc**2)
  call set(sub(boygr,':',':',1), 0.0)

  call grid_bind(l0, 7)
  kappa_l0 = kappa * l0 

  l=q2lb / q2b

  !call tic(2)
  !l=(z>-0.5d0)* max(l, kappa_l0)+(z<=-0.5d0)* l
  ! where(z%data > -0.5) &
  !      l%data = max(l%data, kappa_l0%data)

  call set(l, max(l, kappa_l0), (z > -0.5d0))

  ! gh=l**2 * boygr /q2b
  ! call set(gh, 0.028d0, gh>0.028d0)
  ! gh = l
  ! gh%data = l%data**2 * boygr%data / q2b%data
  ! where(gh%data > 0.028d0) gh%data = 0.028d0

  gh = l**2 * boygr / q2b
  call set(gh, 0.028d0, gh > 0.028d0)
  
  L(1)  = kappa * L0(1)
  L(kb) = 0.d0
  
!  call set(sub(l,':',':',1), kappa*L0(1))
!  call set(sub(l,':',':',kb), 0.d0)
  call set(sub(gh,':',':',1),  0.d0)
  call set(sub(gh,':',':',kb), 0.d0)      

  ! kn= sef*km*(DZB(axf_u) **2 + DZB(ayf_v)**2)/(dhf**2)&
  !      -shiw*km*boygr + kh*boygr

  kn = sef*km*(DZB(AXF(u)) **2 + DZB(AYF(v))**2)/(dhf**2)&
       -shiw*km*boygr + kh*boygr
  

  !     stf=ones(im,jm,kb);
  dtef=sqrt(q2b)/(b1*l+small)
  !      dtef(:,:,1)=0.e0;   dtef(:,:,kb)=0.e0;
  ! DTEF(1) = 0.d0
  ! DTEF(kb)=0.d0
  
  call set(DTEF(1),  0.d0)
  call set(DTEF(kb), 0.d0)

  
  do k=2,kbm1
  !   call set(GG(k),1.d0/(A(k)+C(k)*(1.d0-EE(k-1))-2.d0*dti2*DTEF(k)-1.d0))
  !   call set(EE(k),A(k)*GG(k))
  !   call set(GG(k),(-2.d0*dti2*KN(k)+C(k)*GG(k-1)-Q2F(k))*GG(k))
      GG(k)=1.d0/(A(k)+C(k)*(1.d0-EE(k-1))-2.d0*dti2*DTEF(k)-1.d0)
      EE(k)=A(k)*GG(k)
      GG(k)=(-2.d0*dti2*KN(k)+C(k)*GG(k-1)-Q2F(k))*GG(k)
  enddo

  
  do k=kbm1,1,-1
  !   call set(Q2F(k), EE(k)*Q2F(k+1)+GG(k))
     Q2F(k)=EE(k)*Q2F(k+1)+GG(k)
  enddo
  
  call grid_bind(q2f, 7)

  !!-----------------------solve q2lf
  !call set(Q2LF(kb), 0.d0) 
  Q2LF(kb)= 0.d0

!  call set(EE(2), 0.D0)
!  call set(GG(2), 0.D0)

  do k=2,kbm1
!     call set(DTEF(k),DTEF(k)*(1.e0+e2*((1.e0/abs(z1(k)-z1(1))&
!          +1.e0/abs(z1(k)-z1(kb)))* &
!          L(k)/( (h+etf)*kappa))**2))
     DTEF(k)=DTEF(k)*(1.e0+e2*((1.e0/abs(z1(k)-z1(1))&
          +1.e0/abs(z1(k)-z1(kb)))* &
          L(k)/(DH_2D*kappa))**2)
  enddo

  do k=2,kbm1
!     call set(GG(k),1.d0/(A(k)+C(k)*(1.d0-EE(k-1))-dti2*DTEF(k)-1.d0))
!     call set(EE(k),A(k)*GG(k))
!     call set(GG(k),(-dti2*KN(k)*L(k)*e1+C(k)*GG(k-1)-Q2LF(k))*GG(k))
     GG(k)=1.d0/(A(k)+C(k)*(1.d0-EE(k-1))-dti2*DTEF(k)-1.d0)
     EE(k)=A(k)*GG(k)
     GG(k)=(-dti2*KN(k)*L(k)*e1+C(k)*GG(k-1)-Q2LF(k))*GG(k)
  enddo

  do k=kbm1,2,-1
!     call set(Q2LF(k), EE(k)*Q2LF(k+1)+GG(k))  
     Q2LF(k)= EE(k)*Q2LF(k+1)+GG(k)   
  enddo
  
  call grid_bind(q2lf, 7)

  ! !filter = (q2f<=small .or. q2lf<=small)
  ! filter = mat_zeros
  ! where(q2f%data<=small .or. q2lf%data<=small) &
  !      filter%data = 1.0_8

  filter = (q2f <= small) .or. (q2lf <= small)
  
  call set(sub(filter, ':',':',1),  0)
  call set(sub(filter, ':',':',kb), 0)

  !q2f(filter.data) = small
  !call set(q2f, small, filter=filter)
  ! where(filter%data == 1.0_8) &
  !      q2f%data = small
  call set(q2f, small, filter)
  
  !q2lf(filter.data) = 0.1 * dt(filter.data) * small
  !call set(q2lf, 0.1d0 * dt * small, filter = filter)
  ! where(filter%data == 1.0_8) &
  !      q2lf%data = 0.1d0 * dt%data * small

  call set(q2lf, 0.1d0 * dt * small, filter)
  
  stf=ones(im,jm,kb)
  sh=a2*(1.d0-6.d0*a1/b1*stf)/(1.d0-(3.d0*a2*b2/stf+18.d0*a1*a2)*gh)

  sm=(a1*(1.d0-3.d0*c1-6.d0*a1/b1*stf) &
       + sh*(18.d0*a1*a1+9.d0*a1*a2)*gh)/(1.d0-(9.d0*a1*a2)*gh)

  !kn=l*sqrt(abs(q2));
  kn = l%data * sqrt(abs(q2%data))
  
  !kq=(kn*0.41d0*sh+kq)*0.5d0
  kq = (kn * 0.41d0 * sh + kq) * 0.5d0
  
  !km=(kn*sm+km)*0.5d0;
  km = (kn * sm + km) * 0.5d0;
  
  !kh=(kn*sh+kh)*0.5d0
  kh = (kn * sh + kh) * 0.5d0;

  call set(sub(km,':',jm,':'), &
       sub(km,':',jmm1,':')*sub(fsm,':',jm,':')) 
  call set(sub(kh,':',jm,':'), &
       sub(kh,':',jmm1,':')*sub(fsm,':',jm,':'))     
  call set(sub(km,':',1,':'),  &
       sub(km,':', 2  ,':')*sub(fsm,':',1,':')) 
  call set(sub(kh,':',1,':'),  &
       sub(kh,':', 2  ,':')*sub(fsm,':',1,':'))  

  call set(sub(km, im,':',':'), &
       sub(km,imm1,':',':')*sub(fsm, im,':',':')) 
  call set(sub(kh, im,':',':'), &
       sub(kh,imm1,':',':')*sub(fsm, im,':',':'))    
  call set(sub(km, 1, ':',':'), &
       sub(km,  2, ':',':')*sub(fsm,  1,':',':')) 
  call set(sub(kh, 1, ':',':'), &
       sub(kh,  2, ':',':')*sub(fsm,  1,':',':')) 

  call bcond6(q2f, q2lf,q2, q2l)

  call smoth_update(q2f,q2lf,q2,q2l,q2b,q2lb)

  !utau2,l0,p,cc,boygr,gh,dtef,stf
end subroutine
