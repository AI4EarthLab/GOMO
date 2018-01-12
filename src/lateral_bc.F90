! function [tbe,tbw,tbn,tbs,sbe,sbw,sbn,sbs,uabe,uabw,vabn,vabs] = lateral_bc(bc_flag,iint,time,iend,tbe,tbw,tbn,tbs,sbe,sbw,sbn,sbs,uabe,uabw,vabn,vabs)
subroutine lateral_bc(iint)
  use variables
  use config
  use openarray
  implicit none
  integer :: ibc
  real(kind=8) :: fnew, fold, ntime, sbeb, sbef, sbnb, sbnf
  real(kind=8) :: sbsb, sbsf, sbwb, sbwf, tbc, tbeb, &
       tbef, tbnb,tbnf,tbsb
  real(kind=8) :: tbsf, tbwb, tbwf, uabeb, uabef, &
       uabwb, uabwf, vabnb, vabnf
  real(kind=8) :: vabsf, vabsb
  integer :: iint
  
  if(bc_flag) then
    tbc = 30
    ibc = int(tbc * 86400.0/dti)
    ntime = time / tbc
    if(iint == 1) then
       !call read_boundary_conditions_pnetcdf(iint,ibc);
    end if

    if(iint == 1 .or. mod(iint, ibc) == 0) then
      tbwb = tbwf; sbwb = sbwf;
      tbeb = tbef; sbeb = sbef;
      tbnb = tbnf; sbnb = sbnf;
      tbsb = tbsf; sbsb = sbsf;
      vabnb = vabnf; vabsb = vabsf;
      uabwb = uabwf; uabeb = uabef;
      if(iint /= iend) then
         !call read_boundary_conditions_pnetcdf(iint+ibc, ibc)
      end if
    end if
    
    fnew=time/tbc-ntime;
    fold=1.e0-fnew;
    tbw=fold*tbwb+fnew*tbwf;
    sbw=fold*sbwb+fnew*sbwf;
    tbe=fold*tbeb+fnew*tbef;
    sbe=fold*sbeb+fnew*sbef;
    uabe=fold*uabeb+fnew*uabef;
    uabw=fold*uabwb+fnew*uabwf;
    tbn=fold*tbnb+fnew*tbnf;
    sbn=fold*sbnb+fnew*sbnf;
    tbs=fold*tbsb+fnew*tbsf;
    sbs=fold*sbsb+fnew*sbsf;
    vabn=fold*vabnb+fnew*vabnf;
    vabs=fold*vabsb+fnew*vabsf;
  end if
end subroutine
