-----------------------------------
Data files grl_b2_xx_w(o)em_yyy.dat
-----------------------------------

Horizontal resolution 20 km.

Age-layer data for GRIS.
All data are positive or zero.
Flag value (for NaNs from original file): -666. 

-----
WARNING: there is not any white space separating data
values!! This can be seen in the formatting of 

subroutine read_ad_data(), in read_m.

OpenAD is extremely sensitive to read formats for IO, 
and this is the best solution as yet.
-----

Interpolation on the EPSG:3413 grid,
further processing and downsampling by Liz Logan.

See https://bitbucket.org/lizcurrylogan/matlab_utils/src/master/
files:
> project_age_data.m
> wrt_projected_age_data.m

To see exactly how this data set was created.

Reference
---------

MacGregor, J. A., M. A. Fahnestock, G. A. Catania, J. D. Paden, 
S. P. Gogineni, S. K. Young, S. C. Rybarski, A. N. Mabrey, 
B. M. Wagman, and M. Morlighem (2015), 
Radiostratigraphy and age structure of the Greenland Ice Sheet, 
J. Geophys. Res. Earth Surf., 120, doi:10.1002/2014JF003215.

