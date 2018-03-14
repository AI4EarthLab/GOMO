#include "macro.h"
!#include "common.h"

subroutine internal_q()
  use openarray
  use config
  use variables
  implicit none
  real*8 a1,b1,a2,b2,c1,e1,e2,sef,cbcnst,surfl,shiw
  type(array) :: utau2,l0,p,cc,boygr,gh, stf
  type(array) :: filter, AZB_kq
  type(array) :: dh, kn, sh, sm, kappa_l0
  integer :: ierr, k, i

  type(node) tmp
  
  dhf = h + etf;

  q2f= (q2b*dhb-dti2*(-DZB(AZF(w*q2)) &
       + DXF(AXB(q2) *AXB(dt)* AZB(u) &
       -AZB(AXB(aam))*AXB(h)*DXB(q2b)*dum) &
       +DYF(AYB(q2)* AYB(dt)* AZB(v) &
       -AZB(AYB(aam))*AYB(h)*DYB(q2b)*dvm)))/dhf

  q2lf= (q2lb*dhb-dti2*(-DZB(AZF(w*q2l)) &
       +DXF(AXB(q2l)*AXB(dt)* AZB(u) &
       -AZB(AXB(aam))*AXB(h)* DXB(q2lb)*dum ) &
       +DYF(AYB(q2l)*AYB(dt)*AZB(v)  &
       -AZB(AYB(aam))*AYB(h)*DYB(q2lb)*dvm)))/dhf

  a1=0.92d0;      b1=16.6d0;  a2=0.74d0;
  b2=10.1d0;      c1=0.08d0
  e1=1.8d0;       e2=1.33d0;  sef=1.d0;
  cbcnst=100.d0;  surfl=2.d5; shiw=0.d0

  ee=mat_zeros; gg=mat_zeros; l0=mat_zeros;

  do i = 2, kbm1
    call set(A(i), dzz1(i-1)*dz1(i))
  end do

  a=-dti2*(AZF(kq)+umol)/(a*dhf*dhf)
  call set(A(1),  0.d0)

  do i = 2, kbm1
     call set(C(i), dzz1(i-1)*dz1(i-1))
  end do
  
  c=-dti2*(AZB(kq)+umol)/(c*dhf*dhf)
  call set(C(kb),0.d0)

  utau2 = sqrt(AXF(wusurf)*AXF(wusurf) +AYF(wvsurf)*AYF(wvsurf))
  
  call set(GG(1),(15.8d0*cbcnst)**(2.d0/3.d0)*utau2)

  l0 = rep(surfl*utau2/grav, 1, 1, kb)

  call set(sub(q2f,':',':',kb), &
       sqrt(AXF(wubot)*AXF(wubot) +AYF(wvbot)*AYF(wvbot))*(16.6d0**(2.d0/3.d0))*sef)

  p=grav*rhoref*(-zz * mat_ones * h)*1.d-4

  cc=1449.10d0+0.00821d0*p+4.55d0*(t+tbias) &
       -0.045d0*(t+tbias)*(t+tbias) +1.34d0*(s+sbias-35.0d0)
  
  cc=cc/sqrt((1.d0-0.01642d0*p/cc)*(1.d0-0.4d0*p/(cc*cc)))

  call set(sub(cc,':',':',kb), 0.d0)

  q2b =abs(q2b);  q2lb=abs(q2lb);

  boygr=-grav*DZB(rho)/h + grav*grav/AZB(cc*cc)

  l=q2lb / q2b

  call grid_bind(l0, 7)
  kappa_l0 = kappa * l0 

  call set(l, max(l, kappa_l0), (z_3d > -0.5d0))

  gh = l*l * boygr / q2b
  
  call set(gh, 0.028d0, gh > 0.028d0)

  call set(L(1),  kappa * L0(1))
  call set(L(kb) , 0.d0)
  
  call set(sub(gh,':',':',1),  0.d0)
  call set(sub(gh,':',':',kb), 0.d0)      

  kn = km*sef*(DZB(AXF(u)) *DZB(AXF(u)) + DZB(AYF(v))*DZB(AYF(v)))/(dhf*dhf)&
       -shiw*km*boygr + kh*boygr
  dtef=sqrt(q2b)/(b1*l+small)
  
  call set(DTEF(1),  0.d0)
  call set(DTEF(kb), 0.d0)
  
  do k=2,kbm1
     call set(GG(k),1.d0/(A(k)+C(k)*(1.d0-EE(k-1))-2.d0*dti2*DTEF(k)-1.d0))
     call set(EE(k),A(k)*GG(k))
     call set(GG(k),(-2.d0*dti2*KN(k)+C(k)*GG(k-1)-Q2F(k))*GG(k))
  enddo

  
  do k=kbm1,1,-1
     call set(Q2F(k), EE(k)*Q2F(k+1)+GG(k))
  enddo
  
  call grid_bind(q2f, 7)

  !!-----------------------solve q2lf
  call set(EE(2), 0.d0)
  call set(GG(2), 0.d0)
  call set(Q2LF(kb), 0.d0) 
  call set(Q2LF(1), 0.d0)

  do k=2,kbm1
     call set(DTEF(k),DTEF(k)*(1.d0+e2*((1.d0/abs(z1(k)-z1(1))+1.d0/abs(z1(k)-z1(kb)))*L(k)/(dhf*kappa))*((1.d0/abs(z1(k)-z1(1))+1.d0/abs(z1(k)-z1(kb)))*L(k)/(dhf*kappa))))
  enddo

  do k=2,kbm1
     call set(GG(k),1.d0/(A(k)+C(k)*(1.d0-EE(k-1))-dti2*DTEF(k)-1.d0))
     call set(EE(k),A(k)*GG(k))
     call set(GG(k),(-dti2*KN(k)*L(k)*e1+C(k)*GG(k-1)-Q2LF(k))*GG(k))
  enddo

  do k=kbm1,2,-1
     call set(Q2LF(k), EE(k)*Q2LF(k+1)+GG(k))  
  enddo
  
  call grid_bind(q2lf, 7)

  filter = (q2f <= small) .or. (q2lf <= small)
  
  call set(sub(filter, ':',':',1),  0)
  call set(sub(filter, ':',':',kb), 0)

  call set(q2f, small, filter)

  dt_3d = rep(dt, 1, 1, kb)
  call set(q2lf, 0.1d0 * dt_3d * small, filter)
  
  sh=a2*(1.d0-6.d0*a1/b1)/(1.d0-(3.d0*a2*b2+18.d0*a1*a2)*gh)

  sm=(a1*(1.d0-3.d0*c1-6.d0*a1/b1) &
       + sh*(18.d0*a1*a1+9.d0*a1*a2)*gh)/(1.d0-(9.d0*a1*a2)*gh)

  kn = l * sqrt(abs(q2))
  kq = (kn * 0.41d0 * sh + kq) * 0.5d0
  km = (kn * sm + km) * 0.5d0;
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


  call bcond6()
  call smoth_update()

end subroutine
