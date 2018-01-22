#include "common.h"

module variables
  use openarray
  
  type(array) :: dt,    tbe

  real(kind=8)    rfe, rfw, rfn, rfs

  type(array) aam2d, advua, advva, adx2d, ady2d
  type(array) art, aru, arv, cbc, cor, d, drx2d
  type(array) dry2d   , dum     , dvm
  type(array) dx      , dy      , east_c  , east_e
  type(array) east_u  , east_v  , e_atmos , egb
  type(array) egf     , el      , elb     , elf
  type(array) et      , etb     , etf     , fluxua
  type(array) fluxva  , fsm     , h       , north_c
  type(array) north_e , north_u , north_v , psi
  type(array) rot     , ssurf   , swrad   , vfluxb
  type(array) tps     , tsurf   , ua      , vfluxf
  type(array) uab     , uaf     , utb     , utf
  type(array) va      , vab     , vaf
  type(array) vtb     , vtf     , wssurf  , wtsurf
  type(array) wubot   , wusurf  , wvbot   , wvsurf
  type(array) wubot1  , wvbot1

  type(array) aam  , advx , advy , a
  type(array) c    , drhox, drhoy, dtef
  type(array) ee   , gg   , kh   , km
  type(array) kq   , l    , q2b  , q2, q2f, q2lf
  type(array) q2lb , q2l  , rho  , rmean
  type(array) sb   , sclim, s    , tb
  type(array) tclim, t    , ub   , uf
  type(array) u    , vb   , vf   , v
  type(array) w

  type(array) ele      , eln      , els   , elw
  type(array) sbe,  sbn     , sbs     , sbw
  type(array) tbn     , tbs     , tbw
  type(array) uabe     , uabw     , ube     , ubw
  type(array) vabn     , vabs     , vbn     , vbs

  type(array) z, zz, dz, dzz

  type(array) xflux, yflux, zflux, curv, curv2d, curv_x, curv_y

  ! type(array) :: dx_3d, dy_3d, cor_3d, h_3d, art_3d, aru_3d, arv_3d
  ! type(array) :: fsm_3d, dum_3d, dvm_3d, dt_3d, d_3d

  ! type(array) z_3d, zz_3d, dz_3d, dzz_3d  
  ! type(array)   :: etb_3d,  egf_3d, egb_3d
  ! type(array)   :: e_atmos_3d, etf_3d, dhf_3d, dhb_3d
  ! type(array)   :: axbdt_3d,aybdt_3d

  type(array) z_3d, dt_3d

  type(array)  :: dhf, dhb
  
  type(array)   :: mat_zeros, mat_ones, axf_u,ayf_v

  type(array)  ::  mat_zeros_im_jm_1 
  type(array)  ::  mat_zeros_im_1_1  
  type(array)  ::  mat_zeros_im_1_kb 
  type(array)  ::  mat_zeros_1_jm_kb 
  type(array)  ::  mat_zeros_1_jm_1  
  type(array)  ::  mat_zeros_1_1_kb  
  type(array)  ::  swrad0

  real(kind=8) :: dti, dte2, dti2, ispi, isp2i
  real(kind=8) :: period, time0, time
  integer :: iprint, iswtch, iend
  type(array) :: tf, sf, rad
  
  real(kind=8), allocatable :: z1(:), zz1(:), dzz1(:), dz1(:)

  type(array) :: dxb_axf_dy, dyb_ayf_dx
contains
  subroutine init_variables()
    use config
    implicit none
    character(len=1000) :: fnc
    
    call LoadConfig()
    
    mat_zeros         = zeros(im, jm, kb, sw=1, dt=2)
    mat_ones          = ones(im, jm, kb, sw=1, dt=2)
    mat_zeros_im_jm_1 = sub(mat_zeros, ':', ':', 1)
    mat_zeros_im_1_1  = sub(mat_zeros, ':', 1, 1)
    mat_zeros_im_1_kb = sub(mat_zeros, ':', 1, ':')
    mat_zeros_1_jm_kb = sub(mat_zeros, 1, ':', ':')
    mat_zeros_1_jm_1  = sub(mat_zeros, 1, ':', 1)
    mat_zeros_1_1_kb  = sub(mat_zeros, 1, 1,   ':')
    
    aam=mat_zeros  ;advx=mat_zeros ;advy=mat_zeros ;a=mat_zeros    ;
    c=mat_zeros    ;drhox=mat_zeros;drhoy=mat_zeros;dtef=mat_zeros ;
    ee=mat_zeros   ;gg=mat_zeros   ;kh=mat_zeros   ;km=mat_zeros   ;
    kq=mat_zeros   ;l=mat_zeros    ;q2b=mat_zeros  ;q2=mat_zeros   ;
    q2lb=mat_zeros ;q2l=mat_zeros  ;rho=mat_zeros  ;rmean=mat_zeros;
    sb=mat_zeros   ;sclim=mat_zeros;s=mat_zeros    ;tb=mat_zeros   ;
    tclim=mat_zeros;t=mat_zeros    ;ub=mat_zeros   ;uf=mat_zeros   ;
    u=mat_zeros    ;vb=mat_zeros   ;vf=mat_zeros   ;v=mat_zeros    ;
    w=mat_zeros    ;zflux=mat_zeros;q2f=mat_zeros  ;q2lf=mat_zeros;
    tf = mat_zeros; sf = mat_zeros;
    
    aam2d=mat_zeros_im_jm_1   ;advua=mat_zeros_im_jm_1   ;
    ady2d=mat_zeros_im_jm_1   ;art=mat_zeros_im_jm_1     ;
    cbc=mat_zeros_im_jm_1     ;
    cor=mat_zeros_im_jm_1     ;
    d=mat_zeros_im_jm_1       ;drx2d=mat_zeros_im_jm_1   ;
    dry2d=mat_zeros_im_jm_1   ;dt=mat_zeros_im_jm_1      ;
    dx=mat_zeros_im_jm_1      ;dy=mat_zeros_im_jm_1      ;
    east_c=mat_zeros_im_jm_1  ;east_e=mat_zeros_im_jm_1  ;
    east_u=mat_zeros_im_jm_1  ;east_v=mat_zeros_im_jm_1  ;
    egf=mat_zeros_im_jm_1     ;el=mat_zeros_im_jm_1      ;
    et=mat_zeros_im_jm_1      ;etb=mat_zeros_im_jm_1     ;
    fluxva=mat_zeros_im_jm_1  ;fsm=mat_zeros_im_jm_1     ;
    north_e=mat_zeros_im_jm_1 ;north_u=mat_zeros_im_jm_1 ;
    rot=mat_zeros_im_jm_1     ;ssurf=mat_zeros_im_jm_1   ;
    tps=mat_zeros_im_jm_1     ;tsurf=mat_zeros_im_jm_1   ;
    uab=mat_zeros_im_jm_1     ;uaf=mat_zeros_im_jm_1     ;
    va=mat_zeros_im_jm_1      ;vab=mat_zeros_im_jm_1     ;
    vtb=mat_zeros_im_jm_1     ;vtf=mat_zeros_im_jm_1     ;
    wubot=mat_zeros_im_jm_1   ;wusurf=mat_zeros_im_jm_1  ;
    wubot1=mat_zeros_im_jm_1  ;wvbot1=mat_zeros_im_jm_1  ;


    advva=mat_zeros_im_jm_1   ;adx2d=mat_zeros_im_jm_1   ;
    aru=mat_zeros_im_jm_1     ;arv=mat_zeros_im_jm_1     ;
    dum=mat_zeros_im_jm_1     ;dvm=mat_zeros_im_jm_1     ;
    e_atmos=mat_zeros_im_jm_1 ;egb=mat_zeros_im_jm_1     ;
    elb=mat_zeros_im_jm_1     ;elf=mat_zeros_im_jm_1     ;
    etf=mat_zeros_im_jm_1     ;fluxua=mat_zeros_im_jm_1  ;
    h=mat_zeros_im_jm_1       ;north_c=mat_zeros_im_jm_1 ;
    north_v=mat_zeros_im_jm_1 ;psi=mat_zeros_im_jm_1     ;    
    swrad=mat_zeros_im_jm_1   ;vfluxb=mat_zeros_im_jm_1  ;
    swrad0 = mat_zeros_im_jm_1
    ua=mat_zeros_im_jm_1      ;vfluxf=mat_zeros_im_jm_1  ;
    utb=mat_zeros_im_jm_1     ;utf=mat_zeros_im_jm_1     ;
    vaf=mat_zeros_im_jm_1     ;                         
    wssurf=mat_zeros_im_jm_1  ;wtsurf=mat_zeros_im_jm_1  ;
    wvbot=mat_zeros_im_jm_1   ;wvsurf=mat_zeros_im_jm_1  ;

    ele=mat_zeros_1_jm_1      ;eln=mat_zeros_im_1_1      ;
    sbe=mat_zeros_1_jm_kb     ;sbn=mat_zeros_im_1_kb     ;
    tbe=mat_zeros_1_jm_kb     ;tbn=mat_zeros_im_1_kb     ;
    uabe=mat_zeros_1_jm_1     ;uabw=mat_zeros_1_jm_1     ;
    vabn=mat_zeros_im_1_1     ;vabs=mat_zeros_im_1_1     ;

    els=mat_zeros_im_1_1      ;elw=mat_zeros_1_jm_1      ;
    sbs=mat_zeros_im_1_kb     ;sbw=mat_zeros_1_jm_kb     ;
    tbs=mat_zeros_im_1_kb     ;tbw=mat_zeros_1_jm_kb     ;
    ube=mat_zeros_1_jm_kb     ;ubw=mat_zeros_1_jm_kb     ;
    vbn=mat_zeros_im_1_kb     ;vbs=mat_zeros_im_1_kb     ;

    d   = mat_zeros_im_jm_1;
    dhf = mat_zeros_im_jm_1;
    dhb = mat_zeros_im_jm_1;
    
    !------------------------------------------------------
    !                   start of 3d
    !------------------------------------------------------

    ! d_3d=mat_zeros            ;dt_3d=mat_zeros           ;
    ! dhf_3d=mat_zeros          ;dhb_3d=mat_zeros          ;
    ! egf_3d=mat_zeros          ;etf_3d=mat_zeros          ;
    ! etb_3d=mat_zeros          ;egb_3d=mat_zeros          ;
    ! e_atmos_3d=mat_zeros      ;cor_3d=mat_zeros          ;
    ! z_3d=mat_zeros            ;zz_3d=mat_zeros           ;
    ! dz_3d=mat_zeros           ;dzz_3d=mat_zeros          ;
    ! dx_3d=mat_zeros           ;dy_3d =mat_zeros          ;
    
    ! art_3d=mat_zeros          ;arv_3d=mat_zeros          ;
    ! aru_3d=mat_zeros          ;fsm_3d=mat_zeros          ;
    ! dum_3d=mat_zeros          ;dvm_3d=mat_zeros          ;
    ! h_3d  =mat_zeros          ;

    !------------------------------------------------------
    !                   end of 3d
    !------------------------------------------------------
    
    rad = mat_zeros           

    curv = mat_zeros
    curv_x = mat_zeros;  curv_y = mat_zeros
    curv2d = mat_zeros_im_jm_1

    !call file2ic()
    fnc = trim(in_path)//trim(problem)//".nc"
    print*, "trim(fnc)=", trim(fnc)
    call read_init(trim(fnc))

    dti=dte*isplit;
    
    dte2=dte*2.d0;
    dti2=dti*2.d0;

    !print*, "days = ", days
    iend=max(int(days*24.0*3600.0/dti+0.5),2);
    
    iprint=int(prtd1*24.0*3600.0/dti+0.5);
    iswtch=int(swtch*24.0*3600.0/dti+0.5);

    ispi=1.0/isplit;
    isp2i=1.0/(2.0*isplit);
    !write(*,*) "===end of init_variables"
  end subroutine

end module
