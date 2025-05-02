% =========================================================================
% surfvel_measures_ant_to_sico
%
% Description:
%   Reads the MEaSUREs surface velocity data for the
%   Antarctic ice sheet, downsamples them on the resolutions used by
%   SICOPOLIS and writes them on files in NetCDF format.
%
% Author: Ralf Greve
% Date:   2025-05-02
% =========================================================================

clear variables
close all

%-------- Parameter settings --------

inpath   = '/uchi/greve/sicopolis/data/ant_surfvel_measures';
% inpath   = '.';
filename = 'antarctica_ice_velocity_450m_v2.nc';

xmin = -3040;   % Domain limits of the SICOPOLIS EPSG:3031 grid [km]
xmax =  3040;
ymin = -3040;
ymax =  3040;

ch_grid = 'EPSG3031';

pi_180     = pi/180.0;
pi_180_inv = 180.0/pi;

A = 6.378137e+03;
%   Semi-major axis of the Earth = 6378137 m (in km)

B = 6.3567523142e+03;
%   Semi-minor axis of the Earth = 6356752.3142 m (in km)

lond0 =     0;   % Central meridian  EPSG:3031 [deg]
latd0 =   -71;   % Standard parallel EPSG:3031 [deg]
lon0  = lond0 *pi_180;
lat0  = latd0 *pi_180;

ch_reference = ['Rignot et al., 2017, ' ...
                'NASA National Snow and Ice Data Center ' ...
                'Distributed Active Archive Center, ' ...
                'doi: 10.5067/D7GK8F5J8M8R'];

%-------- Reading of the MEaSUREs surface velocity data --------

disp(' ')
disp(['Reading data file ' filename ' ...'])

ncid = netcdf.open([inpath '/' filename], 'nowrite');

varid = netcdf.inqVarID(ncid, 'x');
x_ms  = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'y');
y_ms  = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'VX');
vx_ms = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'VY');
vy_ms = netcdf.getVar(ncid, varid);

netcdf.close(ncid);

y_ms  = flipud(y_ms);
vx_ms = fliplr(vx_ms);
vy_ms = fliplr(vy_ms);

imax_ms = length(x_ms)-1;
jmax_ms = length(y_ms)-1;

xmin_ms = double(x_ms(1))   *1.0e-03;   % m -> km
xmax_ms = double(x_ms(end)) *1.0e-03;   % m -> km
ymin_ms = double(y_ms(1))   *1.0e-03;   % m -> km
ymax_ms = double(y_ms(end)) *1.0e-03;   % m -> km

dx_ms   = double(x_ms(2)-x_ms(1)) *1.0e-03;   % m -> km

x_ms    = double(x_ms) *1.0e-03;   % m -> km
y_ms    = double(y_ms) *1.0e-03;   % m -> km

vx_ms   = double(vx_ms);
vy_ms   = double(vy_ms);

vs_ms   = sqrt(vx_ms.^2 + vy_ms.^2);

vs_ms(vs_ms<eps) = NaN;

flag_nan = isnan(vs_ms);

%-------- Definition of SICOPOLIS grid --------

disp(' ')
disp('Defining SICOPOLIS grid ...')

disp(' ')
disp('Horizontal resolution: (1) 64 km, (2) 32 km, (3) 16 km, (4) 8 km,')
disp('                       (5) 40 km, (6) 20 km, (7) 10 km.')
dx_number = input('Enter 1, 2, 3, 4, 5, 6 or 7 > ');

if ( (dx_number < 1) || (dx_number > 7) )
   error('Wrong value!')
end

if (dx_number == 1)
   dx    =  64 ;
   ch_dx = '64';
elseif (dx_number == 2)
   dx    =  32 ;
   ch_dx = '32';
elseif (dx_number == 3)
   dx    =  16 ;
   ch_dx = '16';
elseif (dx_number == 4)
   dx    =   8 ;
   ch_dx = '08';
elseif (dx_number == 5)
   dx    =  40 ;
   ch_dx = '40';
elseif (dx_number == 6)
   dx    =  20 ;
   ch_dx = '20';
elseif (dx_number == 7)
   dx    =  10 ;
   ch_dx = '10';
end

disp(' ')
disp(['Resolution = ' ch_dx ' km ...'])

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

filename = ['SurfVel_Antarctica_MEaSUREs_Grid' ch_grid ...
            '_' ch_dx 'km.nc'];

ncid = netcdf.create(filename, 'CLASSIC_MODEL');

dimid_x = netcdf.defDim(ncid, 'x', imax+1);
dimid_y = netcdf.defDim(ncid, 'y', jmax+1);

varid1 = netcdf.defVar(ncid, 'crs', 'NC_BYTE', []);

netcdf.putAtt(ncid, varid1, 'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid, varid1, 'ellipsoid', 'WGS84');
netcdf.putAtt(ncid, varid1, 'false_easting', 0);
netcdf.putAtt(ncid, varid1, 'false_northing', 0);
netcdf.putAtt(ncid, varid1, 'latitude_of_projection_origin', -90);
netcdf.putAtt(ncid, varid1, 'straight_vertical_longitude_from_pole', lond0);
netcdf.putAtt(ncid, varid1, 'standard_parallel', latd0);

varid2 = netcdf.defVar(ncid, 'x', 'NC_INT', dimid_x);

netcdf.putAtt(ncid, varid2, 'units', 'km');
netcdf.putAtt(ncid, varid2, 'standard_name', 'projection_x_coordinate');
netcdf.putAtt(ncid, varid2, 'long_name', 'x-coordinate');
netcdf.putAtt(ncid, varid2, 'axis', 'x');

varid3 = netcdf.defVar(ncid, 'y', 'NC_INT', dimid_y);

netcdf.putAtt(ncid, varid3, 'units', 'km');
netcdf.putAtt(ncid, varid3, 'standard_name', 'projection_y_coordinate');
netcdf.putAtt(ncid, varid3, 'long_name', 'y-coordinate');
netcdf.putAtt(ncid, varid3, 'axis', 'y');

varid4 = netcdf.defVar(ncid, 'lon', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid4, 'units', 'degrees_E');
netcdf.putAtt(ncid, varid4, 'standard_name', 'longitude');
netcdf.putAtt(ncid, varid4, 'long_name', 'Geographical longitude');

varid5 = netcdf.defVar(ncid, 'lat', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid5, 'units', 'degrees_N');
netcdf.putAtt(ncid, varid5, 'standard_name', 'latitude');
netcdf.putAtt(ncid, varid5, 'long_name', 'Geographical latitude');

varid6 = netcdf.defVar(ncid, 'vs', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid6, 'units', 'm a-1');
netcdf.putAtt(ncid, varid6, 'standard_name', ...
                               'land_ice_surface_velocity');
netcdf.putAtt(ncid, varid6, 'long_name', 'Surface velocity');
netcdf.putAtt(ncid, varid6, 'grid_mapping', 'crs');
netcdf.putAtt(ncid, varid6, 'coordinates', 'lon lat');

netcdf.endDef(ncid);

crs = 0;

netcdf.putVar(ncid, varid1, crs);
netcdf.putVar(ncid, varid2, round(x));
netcdf.putVar(ncid, varid3, round(y));
netcdf.putVar(ncid, varid4, mod(lon*pi_180_inv+180, 360)-180);
netcdf.putVar(ncid, varid5, lat*pi_180_inv);
netcdf.putVar(ncid, varid6, vs);

netcdf.close(ncid);

ncwriteatt(filename, '/', 'Conventions', 'CF-1.4');
ncwriteatt(filename, '/', 'Author'     , ...
                          'Ralf Greve <greve@lowtem.hokudai.ac.jp>');
ncwriteatt(filename, '/', 'Variable'   , ...
                          'Surface velocity (m/a) for Antarctica');
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
