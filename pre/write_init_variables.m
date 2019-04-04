file_name = ['./',problem,num2str(im),'_',num2str(jm),'_',num2str(kb),'.nc'];
ncid      =    netcdf.create(file_name,'64BIT_OFFSET');
imx       =    netcdf.defDim(ncid,'im',im);     jmy       =    netcdf.defDim(ncid,'jm',jm);
kbz       =    netcdf.defDim(ncid,'kb',kb);     dconst    =    netcdf.defDim(ncid,'dc',1);
uvel=reshape(uvel,1,1,1);       vvel=reshape(vvel,1,1,1);
z=reshape(z,1,1,kb);            zz=reshape(zz,1,1,kb);              dz=reshape(dz,1,1,kb);          dzz=reshape(dzz,1,1,kb);
dx=reshape(dx,im,jm,1);         dy=reshape(dy,im,jm,1);             cor=reshape(cor,im,jm,1);   
east_e=reshape(east_e,im,jm,1); north_e=reshape(north_e,im,jm,1);     
east_c=reshape(east_c,im,jm,1); north_c=reshape(north_c,im,jm,1);   
east_u=reshape(east_u,im,jm,1); north_u=reshape(north_u,im,jm,1);
east_v=reshape(east_v,im,jm,1); north_v=reshape(north_v,im,jm,1);
h=reshape(h,im,jm,1);           fsm=reshape(fsm,im,jm,1);           dum=reshape(dum,im,jm,1);       dvm=reshape(dvm,im,jm,1);
rot=reshape(rot,im,jm,1);       wusurf=reshape(wusurf,im,jm,1);     wvsurf=reshape(wvsurf,im,jm,1); wtsurf=reshape(wtsurf,im,jm,1);
wtsurf=reshape(wtsurf,im,jm,1); wssurf=reshape(wssurf,im,jm,1);     e_atmos=reshape(e_atmos,im,jm,1);   
swrad=reshape(swrad,im,jm,1);   vflux=reshape(vflux,im,jm,1);       
art=reshape(art,im,jm,1);       aru=reshape(aru,im,jm,1);       arv=reshape(arv,im,jm,1);      

uvelid     =    netcdf.defVar(ncid,'uvel','double',[dconst dconst dconst]);
vvelid     =    netcdf.defVar(ncid,'vvel','double',[dconst dconst dconst]);
zid       =    netcdf.defVar(ncid,'z','double',[dconst dconst kbz]);
zzid      =    netcdf.defVar(ncid,'zz','double',[dconst dconst kbz]);
dxid      =    netcdf.defVar(ncid,'dx','double',[imx jmy dconst]);
dyid      =    netcdf.defVar(ncid,'dy','double',[imx jmy dconst]);
corid     =    netcdf.defVar(ncid,'cor','double',[imx jmy dconst]);
east_eid  =    netcdf.defVar(ncid,'east_e','double',[imx jmy dconst]);
north_eid =    netcdf.defVar(ncid,'north_e','double',[imx jmy dconst]);
east_cid  =    netcdf.defVar(ncid,'east_c','double',[imx jmy dconst]);
north_cid =    netcdf.defVar(ncid,'north_c','double',[imx jmy dconst]);
east_uid  =    netcdf.defVar(ncid,'east_u','double',[imx jmy dconst]);
north_uid =    netcdf.defVar(ncid,'north_u','double',[imx jmy dconst]);
east_vid  =    netcdf.defVar(ncid,'east_v','double',[imx jmy dconst]);
north_vid =    netcdf.defVar(ncid,'north_v','double',[imx jmy dconst]);
hid       =    netcdf.defVar(ncid,'h','double',[imx jmy dconst]);
fsmid     =    netcdf.defVar(ncid,'fsm','double',[imx jmy dconst]);
dumid     =    netcdf.defVar(ncid,'dum','double',[imx jmy dconst]);
dvmid     =    netcdf.defVar(ncid,'dvm','double',[imx jmy dconst]);
rotid     =    netcdf.defVar(ncid,'rot','double',[imx jmy dconst]);
wusurfid  =    netcdf.defVar(ncid,'wusurf','double',[imx jmy dconst]);
wvsurfid  =    netcdf.defVar(ncid,'wvsurf','double',[imx jmy dconst]);
wtsurfid  =    netcdf.defVar(ncid,'wtsurf','double',[imx jmy dconst]);
wssurfid  =    netcdf.defVar(ncid,'wssurf','double',[imx jmy dconst]);
e_atmosid =    netcdf.defVar(ncid,'e_atmos','double',[imx jmy dconst]);
swradid   =    netcdf.defVar(ncid,'swrad','double',[imx jmy dconst]);
vfluxid  =    netcdf.defVar(ncid,'vfluxf','double',[imx jmy dconst]);
tbid     =    netcdf.defVar(ncid,'tb','double',[imx jmy kbz]);
sbid     =    netcdf.defVar(ncid,'sb','double',[imx jmy kbz]);
tclimid  =    netcdf.defVar(ncid,'tclim','double',[imx jmy kbz]);
sclimid  =    netcdf.defVar(ncid,'sclim','double',[imx jmy kbz]);
dzid      =    netcdf.defVar(ncid,'dz','double',[dconst dconst kbz]);
dzzid     =    netcdf.defVar(ncid,'dzz','double',[dconst dconst kbz]);
artid     =    netcdf.defVar(ncid,'art','double',[imx jmy dconst]);
aruid     =    netcdf.defVar(ncid,'aru','double',[imx jmy dconst]);
arvid     =    netcdf.defVar(ncid,'arv','double',[imx jmy dconst]);
ubid     =    netcdf.defVar(ncid,'ub','double',[imx jmy kbz]);
vbid     =    netcdf.defVar(ncid,'vb','double',[imx jmy kbz]);

uab=reshape(uab,im,jm,1);   vab=reshape(vab,im,jm,1);   elb=reshape(elb,im,jm,1);   etb=reshape(etb,im,jm,1);
dt=reshape(dt,im,jm,1);     rfe=reshape(rfe,1,1,1);     rfw=reshape(rfw,1,1,1);     rfn=reshape(rfn,1,1,1);
rfs=reshape(rfs,1,1,1);     uabw=reshape(uabw,1,jm,1);  uabe=reshape(uabe,1,jm,1);  vabs=reshape(vabs,im,1,1);
vabn=reshape(vabn,im,1,1);  elw=reshape(elw,1,jm,1);    ele=reshape(ele,1,jm,1);    els=reshape(els,im,1,1);
eln=reshape(eln,im,1,1);    ssurf=reshape(ssurf,im,jm,1);                           tsurf=reshape(tsurf,im,jm,1);   
tbw=reshape(tbw,1,jm,kb);   tbe=reshape(tbe,1,jm,kb);   tbs=reshape(tbs,im,1,kb);   tbn=reshape(tbn,im,1,kb);
sbw=reshape(sbw,1,jm,kb);   sbe=reshape(sbe,1,jm,kb);   sbs=reshape(sbs,im,1,kb);   sbn=reshape(tbn,im,1,kb);

uabid     =    netcdf.defVar(ncid,'uab','double',[imx jmy dconst]);
vabid     =    netcdf.defVar(ncid,'vab','double',[imx jmy dconst]);
elbid     =    netcdf.defVar(ncid,'elb','double',[imx jmy dconst]);
etbid     =    netcdf.defVar(ncid,'etb','double',[imx jmy dconst]);
dtid      =    netcdf.defVar(ncid,'dt','double',[imx jmy dconst]);
rfeid     =    netcdf.defVar(ncid,'rfe','double',[dconst dconst dconst]);
rfwid     =    netcdf.defVar(ncid,'rfw','double',[dconst dconst dconst]);
rfnid     =    netcdf.defVar(ncid,'rfn','double',[dconst dconst dconst]);
rfsid     =    netcdf.defVar(ncid,'rfs','double',[dconst dconst dconst]);
uabwid     =    netcdf.defVar(ncid,'uabw','double',[dconst jmy dconst]);
uabeid     =    netcdf.defVar(ncid,'uabe','double',[dconst jmy dconst]);
vabsid     =    netcdf.defVar(ncid,'vabs','double',[imx dconst dconst]);
vabnid     =    netcdf.defVar(ncid,'vabn','double',[imx dconst dconst]);
eleid     =    netcdf.defVar(ncid,'ele','double',[dconst jmy dconst]);
elwid     =    netcdf.defVar(ncid,'elw','double',[dconst jmy dconst]);
elsid     =    netcdf.defVar(ncid,'els','double',[imx dconst dconst]);
elnid     =    netcdf.defVar(ncid,'eln','double',[imx dconst dconst]);
ssurfid    =    netcdf.defVar(ncid,'ssurf','double',[imx jmy dconst]);
tsurfid    =    netcdf.defVar(ncid,'tsurf','double',[imx jmy dconst]);
tbeid     =    netcdf.defVar(ncid,'tbe','double',[dconst jmy kbz]);
sbeid     =    netcdf.defVar(ncid,'sbe','double',[dconst jmy kbz]);
sbwid     =    netcdf.defVar(ncid,'sbw','double',[dconst jmy kbz]);
tbwid     =    netcdf.defVar(ncid,'tbw','double',[dconst jmy kbz]);
tbnid     =    netcdf.defVar(ncid,'tbn','double',[imx dconst kbz]);
tbsid     =    netcdf.defVar(ncid,'tbs','double',[imx dconst kbz]);
sbnid     =    netcdf.defVar(ncid,'sbn','double',[imx dconst kbz]);
sbsid     =    netcdf.defVar(ncid,'sbs','double',[imx dconst kbz]);


netcdf.endDef(ncid)
netcdf.putVar(ncid,zid,z);              netcdf.putVar(ncid,zzid,zz);
netcdf.putVar(ncid,dxid,dx);            netcdf.putVar(ncid,dyid,dy);
netcdf.putVar(ncid,corid,cor);          netcdf.putVar(ncid,east_eid,east_e);
netcdf.putVar(ncid,north_eid,north_e);  netcdf.putVar(ncid,east_cid,east_c);
netcdf.putVar(ncid,north_cid,north_c);  netcdf.putVar(ncid,east_uid,east_u);
netcdf.putVar(ncid,north_uid,north_u);  netcdf.putVar(ncid,east_vid,east_v);
netcdf.putVar(ncid,north_vid,north_v);  netcdf.putVar(ncid,hid,h);
netcdf.putVar(ncid,fsmid,fsm);          netcdf.putVar(ncid,rotid,rot);
netcdf.putVar(ncid,dumid,dum);          netcdf.putVar(ncid,dvmid,dvm);
netcdf.putVar(ncid,wusurfid,wusurf);      netcdf.putVar(ncid,wvsurfid,wvsurf);
netcdf.putVar(ncid,wtsurfid,wtsurf);      netcdf.putVar(ncid,wssurfid,wssurf);
netcdf.putVar(ncid,uvelid,uvel);          netcdf.putVar(ncid,vvelid,vvel);
netcdf.putVar(ncid,tbid,tb);              netcdf.putVar(ncid,sbid,sb);
netcdf.putVar(ncid,tclimid,tclim);        netcdf.putVar(ncid,sclimid,sclim);
netcdf.putVar(ncid,e_atmosid,e_atmos);    netcdf.putVar(ncid,swradid,swrad);
netcdf.putVar(ncid,vfluxid,vflux);         netcdf.putVar(ncid,dzid,dz);
netcdf.putVar(ncid,dzzid,dzz);           netcdf.putVar(ncid,artid,art);
netcdf.putVar(ncid,aruid,aru);           netcdf.putVar(ncid,arvid,arv);
netcdf.putVar(ncid,ubid,ub);             netcdf.putVar(ncid,vbid,vb);
netcdf.putVar(ncid,uabid,uab);           netcdf.putVar(ncid,vabid,vab);
netcdf.putVar(ncid,elbid,elb);           netcdf.putVar(ncid,etbid,etb);  
netcdf.putVar(ncid,dtid,dt);             netcdf.putVar(ncid,rfeid,rfe);
netcdf.putVar(ncid,rfwid,rfw);           netcdf.putVar(ncid,rfnid,rfn);
netcdf.putVar(ncid,rfsid,rfs);           netcdf.putVar(ncid,uabwid,uabw); 
netcdf.putVar(ncid,uabeid,uabe);         netcdf.putVar(ncid,vabsid,vabs);
netcdf.putVar(ncid,vabnid,vabn);         netcdf.putVar(ncid,eleid,ele);
netcdf.putVar(ncid,elwid,elw);           netcdf.putVar(ncid,elsid,els);
netcdf.putVar(ncid,elnid,eln);           netcdf.putVar(ncid,ssurfid,ssurf);
netcdf.putVar(ncid,tsurfid,tsurf);       netcdf.putVar(ncid,tbeid,tbe);
netcdf.putVar(ncid,sbeid,sbe);           netcdf.putVar(ncid,sbwid,sbw);
netcdf.putVar(ncid,tbwid,tbw);           netcdf.putVar(ncid,tbnid,tbn);
netcdf.putVar(ncid,tbsid,tbs);           netcdf.putVar(ncid,sbnid,sbn);
netcdf.putVar(ncid,sbsid,sbs);

netcdf.reDef(ncid);
netcdf.putAtt(ncid,zid  ,'ocean_sigma_coordinate_sigma_of_cell_face','dimensionless');
netcdf.putAtt(ncid,zzid ,'ocean_sigma_coordinate_sigma_of_cell_centre','dimensionless');
netcdf.putAtt(ncid,dzid  ,'grid_increment_in_z','dimensionless');
netcdf.putAtt(ncid,dzzid ,'grid_increment_in_zz','dimensionless');

netcdf.putAtt(ncid,dxid ,'grid_increment_in_x','metre');
netcdf.putAtt(ncid,dyid ,'grid_increment_in_y','metre');
netcdf.putAtt(ncid,corid ,'Coriolis_parameter','rad/s');
netcdf.putAtt(ncid,east_eid ,'easting_of_elevation_points','meter');
netcdf.putAtt(ncid,north_eid ,'northing_of_elevation_points','meter');
netcdf.putAtt(ncid,east_cid ,'easting_of_cell_corners','meter');
netcdf.putAtt(ncid,north_cid ,'northing_of_cell_corners','meter');
netcdf.putAtt(ncid,east_uid ,'easting_of_u_point','meter');
netcdf.putAtt(ncid,north_uid ,'northing_of_u_point','meter');
netcdf.putAtt(ncid,east_vid ,'easting_of_v_point','meter');
netcdf.putAtt(ncid,north_vid ,'northing_of_v_point','meter');

netcdf.putAtt(ncid,hid ,'bottom_topography','meter');
netcdf.putAtt(ncid,rotid ,'Rotation_angle_of_x-axis_to_the_east','degree');

netcdf.putAtt(ncid,fsmid ,'free_surface_mask','dimensionless');
netcdf.putAtt(ncid,dumid ,'u-velocity_mask','dimensionless');
netcdf.putAtt(ncid,dvmid ,'v-velocity_mask','dimensionless');

netcdf.putAtt(ncid,wusurfid ,'x-momentum_flux_at_the_surface','metre^2/sec^2');
netcdf.putAtt(ncid,wvsurfid ,'y-momentum_flux_at_the_surface','metre^2/sec^2');
netcdf.putAtt(ncid,wtsurfid ,'temperature_flux_at_the_surface','deg m/s');
netcdf.putAtt(ncid,wssurfid ,'salinity_flux_at_the_surface','psu m/s');

netcdf.putAtt(ncid,uvelid ,'initial_u_velocity','m/s');
netcdf.putAtt(ncid,vvelid ,'initial_v_velocity','m/s');

netcdf.putAtt(ncid,tbid ,'potential_temperature_at_-dt','K');
netcdf.putAtt(ncid,sbid ,'salinity_at_-dt','PSS');
netcdf.putAtt(ncid,tclimid ,'climatology_potential_temperature_at_-dt','K');
netcdf.putAtt(ncid,sclimid ,'climatology_salinity_at_-dt','PSS');

netcdf.putAtt(ncid,e_atmosid ,'depth_changed_by_atmospheric_pressure','metre');
netcdf.putAtt(ncid,swradid,'short_wave_radiation_incident_on_the_ocean_surface','metre/s K');
netcdf.putAtt(ncid,vfluxid,'volumn_flux','metre^3/sec^2');
netcdf.putAtt(ncid,artid,'cell_areas_centered_on_the_variables_t','metre^2');
netcdf.putAtt(ncid,aruid,'cell_areas_centered_on_the_variables_u','metre^2');
netcdf.putAtt(ncid,arvid,'cell_areas_centered_on_the_variables_v','metre^2');

netcdf.putAtt(ncid,ubid,'x-velocity_at_time_-dt','metre/sec');
netcdf.putAtt(ncid,vbid,'y-velocity_at_time_-dt','metre/sec');
netcdf.putAtt(ncid,uabid,'vertical_mean_of_u_at_time_-dt','metre/sec');
netcdf.putAtt(ncid,vabid,'vertical_mean_of_v_at_time_-dt','metre/sec');

netcdf.putAtt(ncid,elbid,'surface_elevation_in_external_mode_at_-dt','metre');
netcdf.putAtt(ncid,etbid,'surface_elevation_in_internal_mode_at_-dt','metre');
netcdf.putAtt(ncid,dtid ,'column_depth_in_internal_mode','metre');

netcdf.putAtt(ncid,rfeid,'Reflection_factor_at_eastern_boundary','dimensionless');
netcdf.putAtt(ncid,rfwid,'Reflection_factor_at_wetern_boundary','dimensionless');
netcdf.putAtt(ncid,rfnid,'Reflection_factor_at_northern_boundary','dimensionless');
netcdf.putAtt(ncid,rfsid,'Reflection_factor_at_southern_boundary','dimensionless');

netcdf.putAtt(ncid,uabeid,'vertical_mean_of_u_at_time_-dt_at_eastern_boundary','metre/sec');
netcdf.putAtt(ncid,uabwid,'vertical_mean_of_u_at_time_-dt_at_western_boundary','metre/sec');
netcdf.putAtt(ncid,vabnid,'vertical_mean_of_v_at_time_-dt_at_northern_boundary','metre/sec');
netcdf.putAtt(ncid,vabsid,'vertical_mean_of_v_at_time_-dt_at_eastern_boundary','metre/sec');

netcdf.putAtt(ncid,eleid,'surface_elevation_in_external_mode_at_eastern_boundary','metre');
netcdf.putAtt(ncid,elwid,'surface_elevation_in_external_mode_at_western_boundary','metre');
netcdf.putAtt(ncid,elnid,'surface_elevation_in_external_mode_at_northern_boundary','metre');
netcdf.putAtt(ncid,elsid,'surface_elevation_in_external_mode_at_sorthern_boundary','metre');

netcdf.putAtt(ncid,ssurfid,'salinity_at_surface','PSS');
netcdf.putAtt(ncid,tsurfid,'potential_temperature','K');

netcdf.putAtt(ncid,tbeid,'potential_temperature_at_time_-dt_at_eastern_boundary','K');
netcdf.putAtt(ncid,tbwid,'potential_temperature_at_time_-dt_at_western_boundary','K');
netcdf.putAtt(ncid,tbnid,'potential_temperature_at_time_-dt_at_northern_boundary','K');
netcdf.putAtt(ncid,tbsid,'potential_temperature_at_time_-dt_at_northern_boundary','K');

netcdf.putAtt(ncid,sbeid,'salinity_at_time_-dt_at_eastern_boundary','PSS');
netcdf.putAtt(ncid,sbwid,'salinity_at_time_-dt_at_western_boundary','PSS');
netcdf.putAtt(ncid,sbnid,'salinity_at_time_-dt_at_northern_boundary','PSS');
netcdf.putAtt(ncid,sbsid,'salinity_at_time_-dt_at_sorthern_boundary','PSS');

netcdf.close(ncid);

disp([file_name,' has been written successully!']);
delete(['seamount_basic',num2str(im),'_',num2str(jm),'_',num2str(kb),'.nc']);
disp('Dimension are ');
im,jm,kb

