% =========================================================================
% surfvel_measures_grl_to_sico
%
% Description:
%   Reads the MEaSUREs surface velocity data for the
%   Greenland ice sheet, downsamples them on the resolutions used by
%   SICOPOLIS and writes them on files in NetCDF format.
%
% Author: Ralf Greve
% Date:   2025-05-02
% =========================================================================

clear variables
close all

%-------- Parameter settings --------

inpath   = '/uchi/greve/sicopolis/data/grl_surfvel_measures';
filename = 'greenland_vel_mosaic250_v1.nc';

xmin =  -720e3;   % Domain limits of the SICOPOLIS EPSG:3413 grid [m]
xmax =   960e3;   % for the GrIS
ymin = -3450e3;
ymax =  -570e3;

ch_grid = 'EPSG3413';

A = 6378137.0;
%   WGS84 semi-major axis [m]

F_INV = 298.257223563;
%   WGS84 inverse flattening

B = A*(1.0-1.0/F_INV);
%   WGS84 semi-minor axis [m]

lond0 = -45;   % Central meridian  EPSG:3413 [deg]
latd0 = +70;   % Standard parallel EPSG:3413 [deg]
lon0  = deg2rad(lond0);
lat0  = deg2rad(latd0);

ch_reference = ['Joughin et al., 2017, ' ...
                'J. Glaciol. 64, 1-11, doi: 10.1017/jog.2017.73'];

%-------- Reading of the MEaSUREs surface velocity data --------

disp(' ')
disp(['Reading data file ' filename ' ...'])

ncid = netcdf.open([inpath '/' filename], 'nowrite');

varid = netcdf.inqVarID(ncid, 'x');
x_ms  = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'y');
y_ms  = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'z');
vs_ms = netcdf.getVar(ncid, varid);

netcdf.close(ncid);

imax_ms = length(x_ms)-1;
jmax_ms = length(y_ms)-1;

xmin_ms = double(x_ms(1));
xmax_ms = double(x_ms(end));
ymin_ms = double(y_ms(1));
ymax_ms = double(y_ms(end));

dx_ms   = double(x_ms(2)-x_ms(1));

x_ms    = double(x_ms);
y_ms    = double(y_ms);

vs_ms   = double(vs_ms);

flag_nan = isnan(vs_ms);

%-------- Definition of SICOPOLIS grid --------

disp(' ')
disp('Defining SICOPOLIS grid ...')

disp(' ')

disp('Horizontal resolution: (1) 40 km, (2) 20 km, (3) 16 km, (4) 10 km,')
disp('                       (5)  8 km, (6)  5 km, (7)  4 km, (8)  2 km.')
dx_number = input('Enter [1-8] > ');

if ( (dx_number < 1) || (dx_number > 8) )
   error('Wrong value!')
end

if (dx_number == 1)
   dx    =  40e3;
   ch_dx = '40';
elseif (dx_number == 2)
   dx    =  20e3;
   ch_dx = '20';
elseif (dx_number == 3)
   dx    =  16e3;
   ch_dx = '16';
elseif (dx_number == 4)
   dx    =  10e3;
   ch_dx = '10';
elseif (dx_number == 5)
   dx    =   8e3;
   ch_dx = '08';
elseif (dx_number == 6)
   dx    =   5e3;
   ch_dx = '05';
elseif (dx_number == 7)
   dx    =   4e3;
   ch_dx = '04';
elseif (dx_number == 8)
   dx    =   2e3;
   ch_dx = '02';
end

disp(' ')
disp(['Resolution = ' ch_dx ' km'])

imax = round((xmax-xmin)/dx);
jmax = round((ymax-ymin)/dx);

x = xmin:dx:xmax;
y = ymin:dx:ymax;

%  ------ Corresponding latitude and longitude

lon = zeros(imax+1,jmax+1);
lat = zeros(imax+1,jmax+1);

for i=1:imax+1
for j=1:jmax+1

   [lon(i,j), lat(i,j)] ...
        = stereo_proj_m('inv_ell', {x(i), y(j), A, B, lon0, lat0});

end
end

%-------- Downsampling of MEaSUREs data on SICOPOLIS grid --------

disp(' ')
disp('Downsampling MEaSUREs data (takes a while) ...')

vs             = zeros(imax+1,jmax+1);
sum_weight     = zeros(imax+1,jmax+1);
sum_weight_max = (dx/dx_ms)^2;

for i=1:imax+1
for j=1:jmax+1

    if (x(i) < xmin_ms || x(i) > xmax_ms) ...
       || ...
       (y(j) < ymin_ms || y(j) > ymax_ms)

        vs(i,j) = 0.0;

    else
        
        x1 = x(i) - 0.5*(dx+dx_ms);
        x2 = x(i) + 0.5*(dx+dx_ms);
        y1 = y(j) - 0.5*(dx+dx_ms);
        y2 = y(j) + 0.5*(dx+dx_ms);

        ii1 = floor((x1-xmin_ms)/dx_ms + 1);
        ii2 =  ceil((x2-xmin_ms)/dx_ms + 1);
        jj1 = floor((y1-ymin_ms)/dx_ms + 1);
        jj2 =  ceil((y2-ymin_ms)/dx_ms + 1);
        
        for ii=ii1:ii2
        for jj=jj1:jj2

            if (ii>=1 && ii<=(imax_ms+1)) ...
               && ...
               (jj>=1 && jj<=(jmax_ms+1))

               if ~isnan(vs_ms(ii,jj))

                   weight = frac_area([x(i)-0.5*dx x(i)+0.5*dx], ...
                                      [y(j)-0.5*dx y(j)+0.5*dx], ...
                              [x_ms(ii)-0.5*dx_ms x_ms(ii)+0.5*dx_ms], ...
                              [y_ms(jj)-0.5*dx_ms y_ms(jj)+0.5*dx_ms]);
                   sum_weight(i,j) = sum_weight(i,j) + weight;
                   vs(i,j)         = vs(i,j) + weight*vs_ms(ii,jj);

               end
               
            end
            
        end
        end
        
        if sum_weight(i,j) > 0.5*sum_weight_max
                             % at least half of the SICOPOLIS grid cell
                             % must be covered by MEaSUREs data
            vs(i,j) = vs(i,j) / sum_weight(i,j);
        else
            vs(i,j) = 0.0;
        end
        
    end
        
end
end

%-------- Writing of data on NetCDF file --------

disp(' ')
disp('Writing data on NetCDF file ...')

filename = ['SurfVel_Greenland_MEaSUREs_Grid' ch_grid ...
            '_' ch_dx 'km.nc'];

ncid = netcdf.create(filename, 'CLASSIC_MODEL');

dimid_x = netcdf.defDim(ncid, 'x', imax+1);
dimid_y = netcdf.defDim(ncid, 'y', jmax+1);

varid1 = netcdf.defVar(ncid, 'mapping', 'NC_BYTE', []);

netcdf.putAtt(ncid, varid1, 'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid, varid1, 'reference_ellipsoid_name', 'WGS84');
netcdf.putAtt(ncid, varid1, 'semi_major_axis', A);
netcdf.putAtt(ncid, varid1, 'inverse_flattening', F_INV);
netcdf.putAtt(ncid, varid1, 'latitude_of_projection_origin', 90);
netcdf.putAtt(ncid, varid1, 'standard_parallel', latd0);
netcdf.putAtt(ncid, varid1, 'straight_vertical_longitude_from_pole', lond0);
netcdf.putAtt(ncid, varid1, 'false_easting', 0);
netcdf.putAtt(ncid, varid1, 'false_northing', 0);

varid2 = netcdf.defVar(ncid, 'x', 'NC_INT', dimid_x);

netcdf.putAtt(ncid, varid2, 'units', 'm');
netcdf.putAtt(ncid, varid2, 'standard_name', 'projection_x_coordinate');
netcdf.putAtt(ncid, varid2, 'long_name', 'x-coordinate');
netcdf.putAtt(ncid, varid2, 'axis', 'x');

varid3 = netcdf.defVar(ncid, 'y', 'NC_INT', dimid_y);

netcdf.putAtt(ncid, varid3, 'units', 'm');
netcdf.putAtt(ncid, varid3, 'standard_name', 'projection_y_coordinate');
netcdf.putAtt(ncid, varid3, 'long_name', 'y-coordinate');
netcdf.putAtt(ncid, varid3, 'axis', 'y');

varid4 = netcdf.defVar(ncid, 'lat', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid4, 'units', 'degrees_N');
netcdf.putAtt(ncid, varid4, 'standard_name', 'latitude');
netcdf.putAtt(ncid, varid4, 'long_name', 'Geographical latitude');

varid5 = netcdf.defVar(ncid, 'lon', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid5, 'units', 'degrees_E');
netcdf.putAtt(ncid, varid5, 'standard_name', 'longitude');
netcdf.putAtt(ncid, varid5, 'long_name', 'Geographical longitude');

varid6 = netcdf.defVar(ncid, 'vs', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid6, 'units', 'm a-1');
netcdf.putAtt(ncid, varid6, 'standard_name', ...
                               'land_ice_surface_velocity');
netcdf.putAtt(ncid, varid6, 'long_name', 'Surface velocity');
netcdf.putAtt(ncid, varid6, 'grid_mapping', 'mapping');

netcdf.endDef(ncid);

mapping = 0;

netcdf.putVar(ncid, varid1, mapping);
netcdf.putVar(ncid, varid2, round(x));
netcdf.putVar(ncid, varid3, round(y));
netcdf.putVar(ncid, varid4, rad2deg(lat));
netcdf.putVar(ncid, varid5, mod(rad2deg(lon)+180, 360)-180);
netcdf.putVar(ncid, varid6, vs);

netcdf.close(ncid);

ncwriteatt(filename, '/', 'Conventions', 'CF-1.4');
ncwriteatt(filename, '/', 'Author'     , ...
                          'Ralf Greve <greve@lowtem.hokudai.ac.jp>');
ncwriteatt(filename, '/', 'Variable'   , ...
                          'Surface velocity (m/a) for Greenland');
ncwriteatt(filename, '/', 'Grid', ch_grid);
ncwriteatt(filename, '/', 'Resolution', [ch_dx ' km']);
ncwriteatt(filename, '/', 'Reference' , ch_reference);
ncwriteatt(filename, '/', 'Date' , datestr(now));

%-------- End of script --------

disp(' ')
disp('Done.')
disp(' ')

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
