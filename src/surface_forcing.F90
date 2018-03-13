#include "common.h"

subroutine surface_forcing(iint)
  use variables
  use config
  use openarray
  implicit none
  real(kind=8) :: fnew, fold, theat, twater, twind
  integer :: iwind, iwater, iheat, ntime, iint
  type(array) :: wtsurfb, wtsurff, swradb, &
       swradf, wssurfb, wssurff, wusurfb, wusurff, &
       wvsurfb, wvsurff
  integer :: ierr

  if (get_rank() == 0) print*, "nsbdy = ", nsbdy, 'iint=',iint
  select case (nsbdy)
  case (1,2,3)
     if(iint == 1) then
       vfluxf = 0.0; wtsurf = 0.0
       satm = 0.0; tatm = 0.0
       wtsurf = mat_zeros_im_jm_1;
       wssurf = mat_zeros_im_jm_1;

     else
        call set(sub(w, ':',':',1), vfluxf)
        return
     endif
  case (4)
     if(wind_flag) then
        twind = 30
        iwind=int(twind*86400.e0/dti);

        if(iwind == 1) then
           call read_var(wusurf, 'wusurf', 'wind', iint, iwind)
           call read_var(wvsurf, 'wvsurf', 'wind', iint, iwind)
        end if
        
        if(iint == 1 .or. mod(iint, iwind) == 0) then
           wusurfb = wusurff; wvsurfb = wvsurff
           if(iint /= iend) then
              call read_var(wusurf, 'wusurf', 'wind', &
                   iint+iwind, iwind)
              call read_var(wvsurf, 'wvsurf', 'wind', &
                   iint+iwind, iwind)
           end if
        end if

        ntime=int(time/twind);
        fnew=time/twind-ntime;
        fold=1.0-fnew;
        wusurf=fold*wusurfb+fnew*wusurff;
        wvsurf=fold*wvsurfb+fnew*wvsurff;
     end if

     if(heat_flag) then
        theat = 1
        iheat=int(theat*86400.e0/dti);
      ! %read heat flux data
      ! % read initial heat file
      if (iint == 1) then
         ![wtsurff,swradf]=read_heat_pnetcdf(iint,iheat);
         call read_var(wtsurff, 'wtsurf', 'heat', iint, iheat)
         call read_var(swradf, 'swrad', 'heat', iint, iheat)
      end if
      
      !% read heat file corresponding to next time
      if (iint == 1 .or. mod(iint,iwind) == 0) then
         wtsurfb=wtsurff;
         swradb=swradf;
         if (iint /= iend) then
            ! [wtsurff,swradf]=read_heat_pnetcdf(iint+iheat,iheat);
            call read_var(wtsurff, 'wtsurf', 'heat', iint, iheat)
            call read_var(swradf, 'swrad', 'heat', iint, iheat)
         endif
      endif

      !% linear interpolation in time
      ntime=int(time/theat);
      fnew=time/theat-ntime;
      fold=1.0-fnew;
      wtsurf=fold*wtsurfb+fnew*wtsurff;
      swrad=fold*swradb+fnew*swradf;
   end if

   if (water_flag) then
      twater=1; !%time between water files(days)

      !%number of steps during a water file      
      iwater=int(twater*86400.e0/dti); 
      
      ! %read water flux data
      ! % read initial water file
      if (iint == 1) then
         ![wssurff]=read_water_pnetcdf(iint,iwater);
         call read_var(wssurff, 'wssurf', 'water', iint, iwater)
      end if
      
      ! % read water file corresponding to next time
      if (iint == 1 .or. mod(iint,iwater) == 0) then
         wssurfb=wssurff;
         if (iint /= iend) then
            ![wssurff]=read_heat_pnetcdf(iint+iwater,iwater);
            call read_var(wssurff, 'wssurf', 'water', iint, iwater) 
         end if
      endif
      
      !% linear interpolation in time
      ntime=int(time/twater);
      fnew=time/twater-ntime;
      fold=1.0-fnew;
      wssurf=fold*wssurfb+fnew*wssurff;
   endif
   
  end select

end subroutine
