
#include "common.h"
#define STENCIL_BOX 1
program pom2k
  use openarray
  use config
  use variables
  implicit none
  integer :: ierr
  integer :: myrank, mysize

  real(kind=8) :: vtot, atot, taver, saver, eaver, tsalt
  integer :: i, iexit, iext, iint, imax, jmax,k
  real(kind=8) :: vamax
  integer :: max_step

  type(array) :: fk, fk1, fk2

  !call tic(11)
  call oa_init(MPI_COMM_WORLD, [-1,-1,1])

  call set_stencil(STENCIL_BOX, 1)

  !get commandline option
  call oa_get_option(max_step, "s", -1)

  !call open_debug()

  !!do i=1,2
  !!do
  !!  !call fuck_you()
  !  fk1 = ones(7,5,6,1,2)
  !  fk2 = zeros(7,5,6,1,2)
  !  fk = fk1 + fk1 - fk1
  !  call disp(fk, "fk = ")
  !!  fk = 1.d0 + fk2 + 1.d0
  !!end do
  !call close_debug()
  !return

  !both namelist data and init data will be loaded.
  call init_variables()

  call grid_init('C', dx, dy, dz)
  !write(*,*) "====end of init_grids"

  call init_fields()
  !write(*,*) "====end of init_fields"

  call update_initial()
  !write(*,*) "====end of update_initial"

  call bottom_friction()
  !write(*,*) "====end of bottom_friction"

  !  open(unit=60,file='conservation7200.txt')
  !  write(60,'(2A21,3A12,A27,A10)')'vtot','atot','eaver','taver','saver','tsalt','vamax'

  do iint = 1, 10

!      if(get_rank()==0) print*,"============iint=============",iint
!      call tic('10')


    call get_time(iint)
    print*, "time = ", time
    print*, "ramp = ", ramp 


    call tic("surf_forcing")
    call surface_forcing(iint)
    call toc("surf_forcing")     
    !call disp(w, "surface_forcing w = ")

    call tic("lateral_bc")     
    call lateral_bc(iint)
    call toc("lateral_bc")          



    call tic("lateral_vis")          
    call lateral_viscosity()
    call toc("lateral_vis")
    !call disp(advx, "advx = ")
    !call disp(advy, "advy = ")
    !call disp(drhox, 'drhox = ')
    !call disp(drhoy, 'drhoy = ')
    !call disp(aam, "aam = ")

    call tic("mode_inter")
    call mode_interaction()
    call toc("mode_inter")
    !call disp(adx2d, "adx2d = ")
    !call disp(ady2d, "ady2d = ")
    !call disp(drx2d, "drx2d = ")
    !call disp(dry2d, "dry2d = ")
    !call disp(aam2d, "aam2d = ")
    !call disp(advua, "advua = ");
    !call disp(advva, "advva = ");
    !call disp(egf, "egf = ")
    !call disp(utf, "utf = ")
    !call disp(vtf, "vtf = ")

     !==============================================
     ! Begin external (2-D) mode
     !==============================================
     do iext=1,int(isplit)
     !do iext=1,30

        call tic("external_el")
        call external_el()
        call toc("external_el")        
        !call disp(elf, "elf = ")


        call tic("external_ua")        
        call external_ua(iext)
        call toc("external_ua")        
        !call disp(uaf, "uaf = ")
        !call disp(advua, "advua = ")

        call tic("external_va")
        call external_va(iext)
        call toc("external_va")        
        !call disp(vaf, "vaf = ")
        !call disp(advva, "advva = ")

        call tic("external_update")
        call external_update(iext, vamax, imax, jmax)
        call toc("external_update")        
        !call disp(utf, "utf = ")
        !call disp(vtf, "vtf = ")

     end do

     fk = vamax * ones(1,1,1,2)
     !call disp(fk, "vamax = ")
     print *, "vamax = ", vamax
     print *, "vmaxl = ", vmaxl
     print *, "isplit = ", int(isplit)

     ! =============================================
     ! End of external (2-D) mode
     ! =============================================

     if(vamax <= vmaxl) then 
        if(((iint/=1).or.(time0/=0.e0) ).and. (mode/=2)) then

           print*,'adjust_uv t='
           call tic("adjust_uv")          
           call adjust_uv()
           call toc("adjust_uv")
           !call disp(u, "u = ")
           !call disp(v, "v = ")

           print*,'internal_w t='
           call tic("internal_w")
           call internal_w()
           call toc("internal_w")
           !call disp(w, "w = ")

           print*,'internal_q t='
           call tic("internal_q")
           call internal_q()
           call toc("internal_q")
           !call disp(q2f, "q2f = ")
           !call disp(q2, "q2 = ")
           !call disp(q2b, "q2b = ")
           !call disp(q2lf, "q2lf = ")
           !call disp(q2l, "q2l = ")
           !call disp(q2lb, "q2lb = ")
           !call disp(km, "km = ")
           !call disp(kq, "kq = ")
           !call disp(kh, "kh = ")
           

           print*,'internal_t t='
           call tic("internal_t1")
           call internal_t(tf,t,tb,wtsurf,tsurf,nbct, &
swrad,tclim,tbe,tbw,tbn,tbs)
           call toc("internal_t1")
           !call disp(tf,'tf = ')
           !call disp(t,'t = ')
           !call disp(tb,'tb =')
           
           print*,'internal_t t='
           call tic("internal_t2")
           call internal_t(sf,s,sb,wssurf,ssurf,nbcs, &
swrad0,sclim,sbe,sbw,sbn,sbs)
           call toc("internal_t2")
           !call disp(sf,'sf = ')
           !call disp(s,'s = ')
           !call disp(sb,'sb =')
           
           print*, "dens"           
           call tic("dens")
           call dens(rho, t, s)
           call toc("dens")
           !call disp(rho, "rho = ")
           
           print*,'internal_u t='
           call tic("internal_u")
           call internal_u()
           call toc("internal_u")
           !call disp(uf, "uf = ")
           !call disp(wubot, "wubot = ")

           print*,'internal_v t='
           call tic("internal_v")
           call internal_v()
           call toc("internal_v")
           !call disp(vf,'vf = ')
           !call disp(wvbot,'wvbot = ')

           print*,'internal_ufvf t='
           call tic("adjust_ufvf")
           call adjust_ufvf()
           call toc("adjust_ufvf")
           !call disp(u, "u = ")
           !call disp(v, "v = ")
           !call disp(ub, "ub = ")
           !call disp(vb, "vb = ")
         endif

         print*, 'internal_update'
         call tic("int_update")
         call internal_update()
         call toc("int_update")

      endif

!      ! print*, "print_seciton"
!      call print_section(iint, iext, vamax,imax,jmax, &
!           vtot, atot, taver, saver, eaver, tsalt)
!      !     write(60,'(F21.7,F21.7,F12.7,F12.7,F16.10,F27.7,F10.5)') &
!      !vtot, atot, eaver,taver,saver ,tsalt,vamax
!      call toc('10')

!      if(iint == max_step) goto 500
!      !     call toc(1)
  enddo

500 print*,"printout at end, iint = ", iint
  !  call disp(uf, 'uf = ')
  !  call disp(vf, 'vf = ')
  !  call disp(advx, 'advx=')
  !  call disp(advy, 'advy=')
  !  call disp(drhox, 'drhox=')
  !  call disp(drhoy, 'drhoy=')
  !  call disp(aam, 'aam =')
  !  call disp(u, 'u =')
  !  call disp(v, 'v =')
  !  call disp(dt_3d, 'dt_3d =')
  !  call disp(w, 'w =')
  !  call disp(q2b,'q2b =')
  !  call disp(q2lb,'q2lb =')
  !  call disp(q2, 'q2 =')
  !  call disp(q2l, 'q2l =')
  !  call disp(q2f, 'q2f =')
  !  call disp(q2lf,'q2lf=')
  !  call disp(km, 'km =')
  !  call disp(kh, 'kh =')
  !  call disp(t, 't =')
  !  call disp(s, 's =')
  !  call disp(rho, 'rho =') 
  !  call disp(ua,  'ua = ')
  !  call disp(va,  'va = ')
  !  call disp(el,  'el = ')
  !  call disp(elf, 'elf = ')
  !  call disp(etf, 'etf = ')
  print*, "vamax = ", vamax
  print*, 'imax=', imax
  print*, 'jmax=', jmax

  ! call save(h, file='dt.nc', var='data')
  ! call system('ncmpidump dt.nc | ncgen -o dt.nc')

  call show_timer()

  call oa_finalize()
end program pom2k
