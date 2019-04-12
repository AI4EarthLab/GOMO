
#include "common.h"
#define STENCIL_BOX 1
program gomo
  use openarray
  use config
  use variables
  use mpi
  implicit none
  integer :: ierr
  integer :: myrank, mysize

  real(kind=8) :: vtot, atot, taver, saver, eaver, tsalt
  integer :: i, iexit, iext, iint, imax, jmax,k
  real(kind=8) :: vamax
  integer :: max_step

  call oa_init(MPI_COMM_WORLD, [-1,-1,1])
  call MPI_COMM_RANK (MPI_COMM_WORLD, myrank, ierr)
  call set_stencil(STENCIL_BOX, 1)

  !get commandline option
  !call oa_get_option(max_step, "s", -1)

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

call tic("zzzall")
  do iint = 1, iend 

    call get_time(iint)

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

    call tic("mode_inter")
    call mode_interaction()
    call toc("mode_inter")

     !==============================================
     ! Begin external (2-D) mode
     !==============================================
     do iext=1,int(isplit)

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
        !call disp(vaf, "external_va vaf = ")
        !call disp(advva, "advva = ")

        call tic("external_update")
        call external_update(iext, vamax, imax, jmax)
        call toc("external_update")        
        !call disp(utf, "utf = ")
        !call disp(vtf, "vtf = ")

     end do

     if (get_rank() == 0) print *, "vamax = ", vamax

     ! =============================================
     ! End of external (2-D) mode
     ! =============================================

     if(vamax <= vmaxl) then 
        if(((iint/=1).or.(time0/=0.e0) ).and. (mode/=2)) then

           !print*,'adjust_uv t='
           call tic("adjust_uv")          
           call adjust_uv()
           call toc("adjust_uv")
           !call disp(u, "u = ")
           !call disp(v, "v = ")

           !print*,'internal_w t='
           call tic("internal_w")
           call internal_w()
           call toc("internal_w")
           !call disp(w, "w = ")

           !print*,'internal_q t='
           call tic("internal_q")
           call internal_q()
           call toc("internal_q")
           
           call tic("internal_t1")
           call internal_t(tf,t,tb,wtsurf,tsurf,nbct, &
swrad,tclim,tbe,tbw,tbn,tbs)
           call toc("internal_t1")
           
           !print*,'internal_t t='
           call tic("internal_t2")
           call internal_t(sf,s,sb,wssurf,ssurf,nbcs, &
swrad0,sclim,sbe,sbw,sbn,sbs)
           call toc("internal_t2")
           
           !print*, "dens"           
           call tic("dens")
           call dens(rho, t, s)
           call toc("dens")
           !call disp(rho, "rho = ")
           
           !print*,'internal_u t='
           call tic("internal_u")
           call internal_u()
           call toc("internal_u")

           !print*,'internal_v t='
           call tic("internal_v")
           call internal_v()
           call toc("internal_v")

           !print*,'internal_ufvf t='
           call tic("adjust_ufvf")
           call adjust_ufvf()
           call toc("adjust_ufvf")

         endif

         !print*, 'internal_update'
         call tic("int_update")
         call internal_update()
         call toc("int_update")

      endif

!      ! print*, "print_seciton"
      call print_section(iint, iext, vamax,imax,jmax, &
           vtot, atot, taver, saver, eaver, tsalt)

  enddo
call toc("zzzall")

  ! call save(h, file='dt.nc', var='data')
  ! call system('ncmpidump dt.nc | ncgen -o dt.nc')

  if(myrank .eq. 0) call show_timer()

  call oa_finalize()
end program gomo
