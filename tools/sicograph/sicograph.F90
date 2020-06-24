!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program   : s i c o g r a p h . F 9 0
!
!> @file
!!
!! Processing and graphical display (with GMT) of the output of SICOPOLIS.
!!
!! @section Date
!!
!! 2020-01-06
!!
!! @section Copyright
!!
!! Copyright 2009-2020 Ralf Greve
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Inclusion of specification header --------

#include RUN_SPECS_HEADER

!-------- Inclusion of kind-type and global-variable modules --------

#include "sicograph_types.F90"
#include "sicograph_variables.F90"

!-------- Further settings --------

#define GMT_SCRIPT_SHELL 'bash'
#define GMT_SCRIPT_PATH './gmt_scripts'

!-------------------------------------------------------------------------------
!> Main program:
!! Processing and graphical display (with GMT) of the output of SICOPOLIS.
!<------------------------------------------------------------------------------
program sicograph

use sicograph_types
use sicograph_variables


implicit none
integer(i4b) :: iexit
integer(i4b) :: menue=0
integer(i4b) :: ndata
character(len=256) :: runname
character(len= 16) :: domain
character(len=  4) :: ergnum


#if (defined(ANT))
   domain = 'ant'
#elif (defined(ASF))
   domain = 'asf'
#elif (defined(EMTP2SGE))
   domain = 'emtp2sge'
#elif (defined(GRL))
   domain = 'grl'
#elif (defined(HEINO))
   domain = 'heino'
#elif (defined(NHEM))
   domain = 'nhem'
#elif (defined(NMARS))
   domain = 'nmars'
#elif (defined(SCAND))
   domain = 'scand'
#elif (defined(SMARS))
   domain = 'smars'
#elif (defined(TIBET))
   domain = 'tibet'
#endif

runname = RUNNAME

!-------- Menue --------

do

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plotting of SICOPOLIS results and data:'
   write(6,'(1x,a)') '---------------------------------------'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plan view:'
   write(6,'(1x,a)') '  (1) Ice-surface topography'
   write(6,'(1x,a)') '  (2) Ice thickness'
   write(6,'(1x,a)') '  (3) Ice-base topography'
   write(6,'(1x,a)') '  (5) Horizontal surface velocity'
   write(6,'(1x,a)') '  (6) Basal temperature (relative to pressure melting)'
   write(6,'(1x,a)') '  (7) Surface temperature'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plan view (data):'
   write(6,'(1x,a)') ' (11) Ice-surface topography'
   write(6,'(1x,a)') ' (12) Ice thickness'
   write(6,'(1x,a)') ' (13) Ice-base topography'
   write(6,'(1x,a)') ' (14) Isostatically relaxed lithosphere surface'
   write(6,'(1x,a)') ' (15) Precipitation rate'
   write(6,'(1x,a)') ' (16) LGM anomaly of precipitation rate'
   write(6,'(1x,a)') ' (17) Surface temperature'
   write(6,'(1x,a)') ' (18) LGM anomaly of surface temperature'
   write(6,'(1x,a)') ' (19) Geothermal heat flux'
   write(6,'(1x,a)') ' (20) Horizontal surface velocity'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') ' (22) (Simulated - measured) ice thickness'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Transect:'
   write(6,'(1x,a)') ' (31) along the x axis [not yet implemented!]'
   write(6,'(1x,a)') ' (32) along the y axis [not yet implemented!]'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Other:'
   write(6,'(1x,a)') ' (41) Time series'
#if (defined(HEINO))
   write(6,'(1x,a)') ' (42) Time series for the sediment region'
#endif
#if ( defined(ANT) \
      || defined(EMTP2SGE) \
      || defined(GRL) \
      || defined(HEINO) \
      || defined(NMARS) )
   write(6,'(1x,a)') ' (43) Time series for ice cores / selected points'
#endif
   write(6,'(1x,a)') ' (44) Scatter plots'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') '  (0) End'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' > '
   read (5,*) menue

   if ((menue >= 1).and.(menue <= 30)) then
      call plan_view(domain, runname, menue)
      menue=0
   else if (menue == 41) then
      call time_series(domain, runname)
      menue=0
#if (defined(HEINO))
   else if (menue == 42) then
      call time_series_sed(domain, runname)
      menue=0
#endif
#if ( defined(ANT) \
      || defined(EMTP2SGE) \
      || defined(GRL) \
      || defined(HEINO) \
      || defined(NMARS) )
   else if (menue == 43) then
      call time_series_core(domain, runname)
      menue=0
#endif
   else if (menue == 44) then
      call scatter_plot(domain, runname)
      menue=0
   end if

   if (menue == 0) exit

end do

write(6,'(1x,a)') ' '

end program sicograph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Plan-view plots (topography etc.).
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine plan_view(domain, runname, menue)

use sicograph_types
use sicograph_variables


implicit none
integer(i4b) :: i, j, n
integer(i4b) :: iexit
integer(i4b) :: ios
integer(i4b) :: menue
integer(i1b) :: i_cb_flag, i_cl_flag
integer(i4b), parameter :: n_unit=11
integer(i1b), dimension(0:IMAX,0:JMAX) :: maske_gr_sim
real(sp), dimension(0:IMAX,0:JMAX) :: zs_gr_sim, zb_gr_sim, &
                                      zl_gr_sim, zl0_gr_sim, &
                                      H_gr_sim
real(sp), dimension(0:IMAX,0:JMAX) :: data
real(sp), parameter :: offset_topo=0.1_sp, offset_thick=0.1_sp, &
                       offset_vel=0.1_sp, offset_temp=0.1_sp, &
                       eps=1.0e-05_sp
character(len=256) :: runname
character(len=256) :: filename_with_path
character(len=256) :: para_file
character(len=256) :: x_min, x_max, x_stp, x_lbl, &
                      y_min, y_max, y_stp, y_lbl, &
                      lambda_min, lambda_max, lambda_stp, lambda_lbl, &
                      phi_min, phi_max, phi_stp, phi_lbl, &
                      dummy_lbl, &
                      lon_ctr, lat_ctr, &
                      lon_lol, lat_lol, lon_upr, lat_upr, &
                      lon_anot, lat_anot, lon_grid, lat_grid
character(len= 16) :: domain, &
                      ch_data, ch_data_unit, dx, dx_int, &
                      dlambda, dlambda_int, dphi, dphi_int, &
                      cb_flag, cl_flag
character(len=  4) :: ergnum

interface
  subroutine read_plot_par(n_unit, lbl, val1, val2, val3, val4, val5)
    use sicograph_types
    implicit none
    integer(i4b) :: n_unit
    character(len=256) :: lbl, val1
    character(len=256), optional :: val2, val3, val4, val5
    end subroutine read_plot_par
end interface


!-------- Grid spacing --------

#if (GRID==0 || GRID==1)

if (abs(DX-4.0).lt.eps) then
   dx =   '4.0'; dx_int =  '04'
else if (abs(DX-5.0).lt.eps) then
   dx =   '5.0'; dx_int =  '05'
else if (abs(DX-10.0).lt.eps) then
   dx =  '10.0'; dx_int =  '10'
else if (abs(DX-20.0).lt.eps) then
   dx =  '20.0'; dx_int =  '20'
else if (abs(DX-25.0).lt.eps) then
   dx =  '25.0'; dx_int =  '25'
else if (abs(DX-40.0).lt.eps) then
   dx =  '40.0'; dx_int =  '40'
else if (abs(DX-50.0).lt.eps) then
   dx =  '50.0'; dx_int =  '50'
else if (abs(DX-75.0).lt.eps) then
   dx =  '75.0'; dx_int =  '75'
else if (abs(DX-80.0).lt.eps) then
   dx =  '80.0'; dx_int =  '80'
else if (abs(DX-100.0).lt.eps) then
   dx = '100.0'; dx_int = '100'
else if (abs(DX-150.0).lt.eps) then
   dx = '150.0'; dx_int = '150'
else if (abs(DX-160.0).lt.eps) then
   dx = '160.0'; dx_int = '160'
else if (abs(DX-200.0).lt.eps) then
   dx = '200.0'; dx_int = '200'
else
   stop ' >>> plan_view: Grid spacing cannot be determined!'
end if

#elif (GRID==2)

if ( (abs(DLAMBDA-(1.0_dp/3.0_dp).lt.eps)) &
     .and. (abs(DPHI-(1.0_dp/3.0_dp).lt.eps)) ) then
   dlambda = '20.0'; dphi = '20.0';  dlambda_int = '20'; dphi_int = '20'
else if ( (abs(DLAMBDA-(1.0_dp/6.0_dp).lt.eps)) &
          .and. (abs(DPHI-(1.0_dp/6.0_dp).lt.eps)) ) then
   dlambda = '10.0'; dphi = '10.0';  dlambda_int = '10'; dphi_int = '10'
else if ( (abs(DLAMBDA-(1.0_dp/12.0_dp).lt.eps)) &
          .and. (abs(DPHI-(1.0_dp/12.0_dp).lt.eps)) ) then
   dlambda = '5.0'; dphi = '5.0';  dlambda_int = '05'; dphi_int = '05'
else
   stop ' >>> plan_view: Grid spacing cannot be determined!'
end if

#endif

!-------- Reading of time-slice data --------

if ((menue >= 1).and.(menue <= 10)) then
   call read_erg_nc(runname, ergnum, menue)
end if

!-------- Reading of input data --------

if ((menue >= 11).and.(menue <= 20)) then

#if (GRID==0 || GRID==1)
   call read_data(domain, dx_int, menue)
#elif (GRID==2)
   call read_data(domain, dphi_int, menue)
#endif
   ergnum = 'data'
end if

!-------- Reading of time-slice and input data for difference plots --------

if ((menue >= 21).and.(menue <= 30)) then

   call read_erg_nc(runname, ergnum, menue)

   zs_gr_sim    = zs_gr
   zb_gr_sim    = zb_gr
   zl_gr_sim    = zl_gr
   zl0_gr_sim   = zl0_gr
   maske_gr_sim = maske_gr
   H_gr_sim     = H_gr

#if (GRID==0 || GRID==1)
   call read_data(domain, dx_int, menue)
#elif (GRID==2)
   call read_data(domain, dphi_int, menue)
#endif

end if

!-------- Plotting of colour bar --------

write(6,'(1x,a)') ' '
write(6,'(1x,a)',advance='no') &
        'Plot (1) with or (2) without colour bar? > '
read (5,*) i_cb_flag

if (i_cb_flag == 1) then
   cb_flag = 'true'
else
   cb_flag = 'false'
end if

!-------- Plotting of contour labels --------

write(6,'(1x,a)') ' '
write(6,'(1x,a)',advance='no') &
        'Plot (1) with or (2) without contour labels? > '
! write(6,'(1x,a)') '(Will be overridden if settings are made in *.zzz file.)'
read (5,*) i_cl_flag

if (i_cl_flag == 1) then
   cl_flag = 'true'
else
   cl_flag = 'false'
end if

!-------- Writing of data on temporary file --------

!  ------ Quantity to be plotted

if ((menue==1).or.(menue==11)) then

   data         =  zs_gr
   ch_data      = 'zs'
   ch_data_unit = 'km'

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   data = data - z_sl_gr   ! topography relative to sea level

   do i=0, IMAX
   do j=0, JMAX
      if (maske_gr(i,j) == 2) data(i,j) = min(data(i,j), -offset_topo)
   end do
   end do

#endif

else if ((menue==2).or.(menue==12)) then

   data         =  H_gr
   ch_data      = 'H'
   ch_data_unit = 'km'

   do i=0, IMAX
   do j=0, JMAX
      if ((maske_gr(i,j) == 1).or.(maske_gr(i,j) == 2)) &
         data(i,j) = -offset_thick
   end do
   end do

else if ((menue==3).or.(menue==13)) then

   data         =  zb_gr
   ch_data      = 'zb'
   ch_data_unit = 'km'

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   data = data - z_sl_gr   ! topography relative to sea level

   do i=0, IMAX
   do j=0, JMAX
      if (maske_gr(i,j) == 2) data(i,j) = min(data(i,j), -offset_topo)
   end do
   end do
#endif

else if (menue==14) then

   data         =  zl0_gr
   ch_data      = 'zl0'
   ch_data_unit = 'km'

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   data = data - z_sl_gr   ! topography relative to sea level

   do i=0, IMAX
   do j=0, JMAX
      if (maske_gr(i,j) == 2) data(i,j) = min(data(i,j), -offset_topo)
   end do
   end do

#endif

else if (menue==5) then

   do i=0, IMAX
   do j=0, JMAX
      if ((maske_gr(i,j) == 0).or.(maske_gr(i,j) == 3)) then
         data(i,j) = vh_s_gr(i,j)
      else
         data(i,j) = -offset_vel
      end if
   end do
   end do

   ch_data      = 'vs'

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   ch_data_unit = 'm/a'

#else   /* Martian ice sheet */

   ch_data_unit = 'mm/a'

#endif

else if (menue==6) then

   do i=0, IMAX
   do j=0, JMAX

      data(i,j) = temph_b_gr(i,j)

/* #if ( defined(ANT) \                                    */
/*       || defined(ASF) \                                 */
/*       || defined(GRL) \                                 */
/*       || defined(SCAND) \                               */
/*       || defined(NHEM) \                                */
/*       || defined(TIBET) \				   */
/*       || defined(EMTP2SGE) \                            */
/*       || defined(HEINO) )   terrestrial ice sheet       */

      if ( (maske_gr(i,j) == 1).or.(maske_gr(i,j) == 2) &
                               .or.(n_cts_gr(i,j) /= -1) ) &
         data(i,j) = offset_temp

/* #endif */

   end do
   end do

   ch_data      = 'Tbh'
   ch_data_unit = 'C'   ! 'deg\ C'

else if (menue==7) then

   do i=0, IMAX
   do j=0, JMAX
      data(i,j) = temp_s_gr(i,j)
   end do
   end do

   ch_data      = 'Ts'
   ch_data_unit = 'C'   ! 'deg\ C'

else if (menue==15) then

   data         =  precip_gr
   ch_data      = 'precip'
   ch_data_unit = 'mm\ w.e./a'

else if (menue==16) then

   data         =  precip_lgm_anom_gr
   ch_data      = 'precip_lgm_anom'
   ch_data_unit = '1'

else if (menue==17) then

   data         =  temp_s_gr
   ch_data      = 'temp_s'
   ch_data_unit = 'C'   ! 'deg\ C'

else if (menue==18) then

   data         =  temp_s_lgm_anom_gr
   ch_data      = 'temp_s_lgm_anom'
   ch_data_unit = 'C'   ! 'deg\ C'

else if (menue==19) then

   data         =  qgeo_gr
   ch_data      = 'qgeo'
   ch_data_unit = 'mW/m2'

else if (menue==20) then

   do i=0, IMAX
   do j=0, JMAX
      data(i,j) =  vh_s_gr(i,j)
      if (data(i,j) < eps) data(i,j) = -offset_vel
   end do
   end do

   ch_data      = 'vs'

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   ch_data_unit = 'm/a'

#else   /* Martian ice sheet */

   ch_data_unit = 'mm/a'

#endif

else if (menue==22) then

   data         = H_gr_sim-H_gr
   data         = data * 1000.0   ! km -> m
   ch_data      = 'dH'
   ch_data_unit = 'm'

end if

filename_with_path = trim(GMT_SCRIPT_PATH)//'/tmp/' &
                     //trim(runname)//trim(ergnum)//'_' &
                     //trim(ch_data)//'_tmp.erg'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do i=0, IMAX
do j=0, JMAX
   write(n_unit,'(3(1x,1pe10.3))') xi_gr(i), eta_gr(j), data(i,j)
end do
end do

close(n_unit, status='keep')

!  ------ Ice thickness (for plotting the ice margin)

data = H_gr

do i=0, IMAX
do j=0, JMAX
   if ((maske_gr(i,j) == 1).or.(maske_gr(i,j) == 2)) &
      data(i,j) = -offset_thick
end do
end do

filename_with_path = trim(GMT_SCRIPT_PATH)//'/tmp/' &
                     //trim(runname)//trim(ergnum)//'_' &
                     //trim(ch_data)//'_tmp2.erg'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do i=0, IMAX
do j=0, JMAX
   write(n_unit,'(3(1x,1pe10.3))') xi_gr(i), eta_gr(j), data(i,j)
end do
end do

close(n_unit, status='keep')

!  ------ Surface topography (for plotting the h=0 contour)

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

data = zs_gr

data = data - z_sl_gr   ! topography relative to sea level

do i=0, IMAX
do j=0, JMAX
   if (maske_gr(i,j) == 2) data(i,j) = min(data(i,j), -offset_topo)
end do
end do

#else   /* Martian ice sheet */

data = -100.0   ! dummy

#endif

filename_with_path = trim(GMT_SCRIPT_PATH)//'/tmp/' &
                     //trim(runname)//trim(ergnum)//'_' &
                     //trim(ch_data)//'_tmp3.erg'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do i=0, IMAX
do j=0, JMAX
   write(n_unit,'(3(1x,1pe10.3))') xi_gr(i), eta_gr(j), data(i,j)
end do
end do

close(n_unit, status='keep')

!-------- Reading of plot parameters --------

para_file          = 'plan_view_'//trim(domain)//'.dat'
filename_with_path = 'parameter_files/'//trim(para_file)

open(n_unit, iostat=ios, file=trim(filename_with_path), status='old')

if (ios == 0) then

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' available.'
   write(6,'(1x,a)') 'Parameters listed in this file used.'

#if (GRID==0 || GRID==1)
   call read_plot_par(n_unit, x_lbl, x_min, x_max, x_stp)
   call read_plot_par(n_unit, y_lbl, y_min, y_max, y_stp)
   call read_plot_par(n_unit, dummy_lbl, lon_ctr)
   call read_plot_par(n_unit, dummy_lbl, lat_ctr)
   call read_plot_par(n_unit, dummy_lbl, lon_lol)
   call read_plot_par(n_unit, dummy_lbl, lat_lol)
   call read_plot_par(n_unit, dummy_lbl, lon_upr)
   call read_plot_par(n_unit, dummy_lbl, lat_upr)
   call read_plot_par(n_unit, dummy_lbl, lon_anot)
   call read_plot_par(n_unit, dummy_lbl, lat_anot)
   call read_plot_par(n_unit, dummy_lbl, lon_grid)
   call read_plot_par(n_unit, dummy_lbl, lat_grid)
#elif (GRID==2)
   call read_plot_par(n_unit, lambda_lbl, lambda_min, lambda_max, lambda_stp)
   call read_plot_par(n_unit, phi_lbl, phi_min, phi_max, phi_stp)
#endif

else

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' not available.'
   write(6,'(1x,a)') 'Parameters will be computed automatically.'

#if (GRID==0 || GRID==1)
   x_lbl = 'none'   ! dummy
   x_min = 'none'   ! dummy
   x_max = 'none'   ! dummy
   x_stp = 'none'   ! dummy

   y_lbl = 'none'   ! dummy
   y_min = 'none'   ! dummy
   y_max = 'none'   ! dummy
   y_stp = 'none'   ! dummy

   lon_ctr  = 'none'   ! dummy
   lat_ctr  = 'none'   ! dummy
   lon_lol  = 'none'   ! dummy
   lat_lol  = 'none'   ! dummy
   lon_upr  = 'none'   ! dummy
   lat_upr  = 'none'   ! dummy
   lon_anot = 'none'   ! dummy
   lat_anot = 'none'   ! dummy
   lon_grid = 'none'   ! dummy
   lat_grid = 'none'   ! dummy
#elif (GRID==2)
   lambda_lbl = 'none'   ! dummy
   lambda_min = 'none'   ! dummy
   lambda_max = 'none'   ! dummy
   lambda_stp = 'none'   ! dummy

   phi_lbl = 'none'   ! dummy
   phi_min = 'none'   ! dummy
   phi_max = 'none'   ! dummy
   phi_stp = 'none'   ! dummy
#endif

end if

close(n_unit, status='keep')

!-------- Plotting of data --------

#if (GRID==0 || GRID==1)

   call system(GMT_SCRIPT_SHELL//' '//GMT_SCRIPT_PATH &
         //'/plan_view.gmt' &
         //' '//trim(domain) &
         //' '//trim(runname)//trim(ergnum)//'_'//trim(ch_data) &
         //' '//trim(ch_data) &
         //' '//trim(ch_data_unit) &
         //' '//trim(adjustl(dx)) &
         //' '//trim(adjustl(x_min)) &
         //' '//trim(adjustl(x_max)) &
         //' '//trim(adjustl(x_stp)) &
         //' '//trim(adjustl(y_min)) &
         //' '//trim(adjustl(y_max)) &
         //' '//trim(adjustl(y_stp)) &
         //' '//trim(adjustl(x_lbl)) &
         //' '//trim(adjustl(y_lbl)) &
         //' '//trim(adjustl(lon_ctr)) &
         //' '//trim(adjustl(lat_ctr)) &
         //' '//trim(adjustl(lon_lol)) &
         //' '//trim(adjustl(lat_lol)) &
         //' '//trim(adjustl(lon_upr)) &
         //' '//trim(adjustl(lat_upr)) &
         //' '//trim(adjustl(lon_anot)) &
         //' '//trim(adjustl(lat_anot)) &
         //' '//trim(adjustl(lon_grid)) &
         //' '//trim(adjustl(lat_grid)) &
         //' '//trim(cb_flag) &
         //' '//trim(cl_flag))

#elif (GRID==2)

   call system(GMT_SCRIPT_SHELL//' '//GMT_SCRIPT_PATH &
         //'/plan_view_lonlat.gmt' &
         //' '//trim(domain) &
         //' '//trim(runname)//trim(ergnum)//'_'//trim(ch_data) &
         //' '//trim(ch_data) &
         //' '//trim(ch_data_unit) &
         //' '//trim(adjustl(dlambda)) &
         //' '//trim(adjustl(dphi)) &
         //' '//trim(adjustl(lambda_min)) &
         //' '//trim(adjustl(lambda_max)) &
         //' '//trim(adjustl(lambda_stp)) &
         //' '//trim(adjustl(phi_min)) &
         //' '//trim(adjustl(phi_max)) &
         //' '//trim(adjustl(phi_stp)) &
         //' '//trim(adjustl(lambda_lbl)) &
         //' '//trim(adjustl(phi_lbl)) &
         //' '//trim(cb_flag) &
         //' '//trim(cl_flag))

#endif

end subroutine plan_view

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Plotting of time-series data (volume, area etc.).
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_series(domain, runname)

use sicograph_types
use sicograph_variables


implicit none
integer(i4b) :: n
integer(i4b) :: ndata
integer(i4b) :: iexit
integer(i4b) :: ios
integer(i4b) :: menue=0
integer(i4b), parameter :: n_unit=11
real(sp), dimension(1048576) :: data
character(len=256) :: runname
character(len=256) :: filename_with_path
character(len=256) :: para_file
character(len=256) :: time_min, time_max, time_stp, time_lbl, &
                      data_min, data_max, data_stp, data_lbl
character(len= 16) :: domain, ch_data
logical :: menue_flag=.false.

interface
  subroutine read_plot_par(n_unit, lbl, val1, val2, val3, val4, val5)
    use sicograph_types
    implicit none
    integer(i4b) :: n_unit
    character(len=256) :: lbl, val1
    character(len=256), optional :: val2, val3, val4, val5
    end subroutine read_plot_par
end interface


!-------- Reading of time-series data --------

call read_ser(runname, ndata)

!-------- Menue --------

do

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') '  (1) Temperature-deviation forcing'
   write(6,'(1x,a)') '  (2) Sea-level forcing'
   write(6,'(1x,a)') '  (3) Maximum ice thickness'
   write(6,'(1x,a)') '  (4) Maximum surface elevation'
   write(6,'(1x,a)') '  (5) Total ice volume'
   write(6,'(1x,a)') '  (6) Temperate ice volume'
#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */
   write(6,'(1x,a)') '  (8) Sea-level equivalent'
#else   /* Martian ice sheet */
   write(6,'(1x,a)') '  (8) Maximum basal temperature (relative to pressure melting)'
#endif
   write(6,'(1x,a)') '  (9) Total ice area' 
   write(6,'(1x,a)') ' (10) Area covered by temperate ice'
   write(6,'(1x,a)') ' (13) Maximum thickness of temperate ice'
   write(6,'(1x,a)') ' (14) Maximum surface velocity'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' > '
   read (5,*) menue

   if (menue == 1) then
      data    =  D_Ts
      ch_data = 'D_Ts'
      menue_flag = .true.
   else if (menue == 2) then
      data    =  z_sl
      ch_data = 'z_sl'
      menue_flag = .true.
   else if (menue == 3) then
      data    =  H_max
      ch_data = 'H_max'
      menue_flag = .true.
   else if (menue == 4) then
      data    =  zs_max
      ch_data = 'zs_max'
      menue_flag = .true.
   else if (menue == 5) then
      data    =  V_tot
      ch_data = 'V_tot'
      menue_flag = .true.
   else if (menue == 6) then
      data    =  V_temp
      ch_data = 'V_temp'
      menue_flag = .true.
   else if (menue == 8) then
#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */
      data    =  V_sle
      ch_data = 'V_sle'
#else   /* Martian ice sheet */
      data    =  Tbh_max
      ch_data = 'Tbh_max'
#endif
      menue_flag = .true.
   else if (menue == 9) then
      data    =  A_tot
      ch_data = 'A_tot'
      menue_flag = .true.
   else if (menue == 10) then
      data    =  A_temp
      ch_data = 'A_temp'
      menue_flag = .true.
   else if (menue == 13) then
      data    =  H_t_max
      ch_data = 'H_t_max'
      menue_flag = .true.
   else if (menue == 14) then
      data    =  vs_max
      ch_data = 'vs_max'
      menue_flag = .true.
   end if

   if (menue_flag) exit

end do

!-------- Writing of data on temporary file --------

filename_with_path = GMT_SCRIPT_PATH//'/tmp/' &
                     //trim(runname)//'_'//trim(ch_data)//'_tmp.ser'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do n=1, ndata
   write(n_unit,'(2(2x,1pe10.3))') time(n), data(n)
end do

close(n_unit, status='keep')

!-------- Reading of plot parameters --------

para_file          = 'time_series_'//trim(domain)//'.dat'
filename_with_path = 'parameter_files/'//trim(para_file)

open(n_unit, iostat=ios, file=trim(filename_with_path), status='old')

if (ios == 0) then

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' available.'
   write(6,'(1x,a)') 'Parameters listed in this file used.'

   call read_plot_par(n_unit, time_lbl, time_min, time_max, time_stp)
   do n=1, menue
      call read_plot_par(n_unit, data_lbl, data_min, data_max, data_stp)
   end do

else

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' not available.'
   write(6,'(1x,a)') 'Parameters will be computed automatically.'

   time_lbl = 'none'   ! dummy
   time_min = 'none'   ! dummy
   time_max = 'none'   ! dummy
   time_stp = 'none'   ! dummy

   data_lbl = 'none'   ! dummy
   data_min = 'none'   ! dummy
   data_max = 'none'   ! dummy
   data_stp = 'none'   ! dummy

end if

close(n_unit, status='keep')

!-------- Plotting of data --------

   call system(GMT_SCRIPT_SHELL//' '//GMT_SCRIPT_PATH &
         //'/time_series.gmt' &
         //' '//trim(domain) &
         //' '//trim(runname)//'_'//trim(ch_data) &
         //' '//trim(ch_data) &
         //' '//trim(adjustl(time_min)) &
         //' '//trim(adjustl(time_max)) &
         //' '//trim(adjustl(time_stp)) &
         //' '//trim(adjustl(data_min)) &
         //' '//trim(adjustl(data_max)) &
         //' '//trim(adjustl(data_stp)) &
         //' '//trim(adjustl(time_lbl)) &
         //' '//trim(adjustl(data_lbl)))

end subroutine time_series

#if (defined(HEINO))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Plotting of time-series data for the sediment region of the HEINO set-up.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_series_sed(domain, runname)

use sicograph_types
use sicograph_variables

#if (defined(NAG))
use f90_unix_proc
    ! Required for executing shell commands with "call system('...')"
#endif

implicit none
integer(i4b) :: n
integer(i4b) :: ndata
integer(i4b) :: iexit
integer(i4b) :: ios
integer(i4b) :: menue=0
integer(i4b), parameter :: n_unit=11
real(sp), dimension(1048576) :: data
character(len=256) :: runname
character(len=256) :: filename_with_path
character(len=256) :: para_file
character(len=256) :: time_min, time_max, time_stp, time_lbl, &
                      data_min, data_max, data_stp, data_lbl
character(len= 16) :: domain, ch_data
logical :: menue_flag=.false.

interface
  subroutine read_plot_par(n_unit, lbl, val1, val2, val3, val4, val5)
    use sicograph_types
    implicit none
    integer(i4b) :: n_unit
    character(len=256) :: lbl, val1
    character(len=256), optional :: val2, val3, val4, val5
    end subroutine read_plot_par
end interface


!-------- Reading of time-series data --------

call read_sed(runname, ndata)

!-------- Menue --------

do

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') '  (1) Temperature-deviation forcing'
   write(6,'(1x,a)') '  (2) Sea-level forcing'
   write(6,'(1x,a)') '  (3) Average ice thickness'
   write(6,'(1x,a)') '  (4) Average basal temperature (relative to pressure melting)'
   write(6,'(1x,a)') '  (5) Temperate basal area'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' > '
   read (5,*) menue

   if (menue == 1) then
      data    =  D_Ts
      ch_data = 'D_Ts'
      menue_flag = .true.
   else if (menue == 2) then
      data    =  z_sl
      ch_data = 'z_sl'
      menue_flag = .true.
   else if (menue == 3) then
      data    =  H_ave_sed
      ch_data = 'H_ave'
      menue_flag = .true.
   else if (menue == 4) then
      data    =  Tbh_ave_sed
      ch_data = 'Tbh_ave'
      menue_flag = .true.
   else if (menue == 5) then
      data    =  Atb_sed
      ch_data = 'Atb'
      menue_flag = .true.
   end if

   if (menue_flag) exit

end do

!-------- Writing of data on temporary file --------

filename_with_path = GMT_SCRIPT_PATH//'/tmp/' &
                     //trim(runname)//'_'//trim(ch_data)//'_tmp.ser'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do n=1, ndata
   write(n_unit,'(2(2x,1pe10.3))') time(n), data(n)
end do

close(n_unit, status='keep')

!-------- Reading of plot parameters --------

para_file          = 'time_series_sed_'//trim(domain)//'.dat'
filename_with_path = 'parameter_files/'//trim(para_file)

open(n_unit, iostat=ios, file=trim(filename_with_path), status='old')

if (ios == 0) then

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' available.'
   write(6,'(1x,a)') 'Parameters listed in this file used.'

   call read_plot_par(n_unit, time_lbl, time_min, time_max, time_stp)
   do n=1, menue
      call read_plot_par(n_unit, data_lbl, data_min, data_max, data_stp)
   end do

else

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' not available.'
   write(6,'(1x,a)') 'Parameters will be computed automatically.'

   time_lbl = 'none'   ! dummy
   time_min = 'none'   ! dummy
   time_max = 'none'   ! dummy
   time_stp = 'none'   ! dummy

   data_lbl = 'none'   ! dummy
   data_min = 'none'   ! dummy
   data_max = 'none'   ! dummy
   data_stp = 'none'   ! dummy

end if

close(n_unit, status='keep')

!-------- Plotting of data --------

call system(GMT_SCRIPT_SHELL//' '//GMT_SCRIPT_PATH &
         //'/time_series.gmt' &
         //' '//trim(domain) &
         //' '//trim(runname)//'_'//trim(ch_data) &
         //' '//trim(ch_data) &
         //' '//trim(adjustl(time_min)) &
         //' '//trim(adjustl(time_max)) &
         //' '//trim(adjustl(time_stp)) &
         //' '//trim(adjustl(data_min)) &
         //' '//trim(adjustl(data_max)) &
         //' '//trim(adjustl(data_stp)) &
         //' '//trim(adjustl(time_lbl)) &
         //' '//trim(adjustl(data_lbl)))

end subroutine time_series_sed

#endif

#if ( defined(ANT) \
      || defined(EMTP2SGE) \
      || defined(GRL) \
      || defined(HEINO) \
      || defined(NMARS) )

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Plotting of time-series data for ice cores or selected points.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine time_series_core(domain, runname)

use sicograph_types
use sicograph_variables

#if (defined(NAG))
use f90_unix_proc
    ! Required for executing shell commands with "call system('...')"
#endif

implicit none
integer(i4b) :: n
integer(i4b) :: ndata
integer(i4b) :: iexit
integer(i4b) :: ios
integer(i4b) :: menue=0, n_core=0
integer(i4b), parameter :: n_unit=11
real(sp), dimension(1048576) :: data
character(len=256) :: runname
character(len=256) :: filename_with_path
character(len=256) :: para_file
character(len=256) :: time_min, time_max, time_stp, time_lbl, &
                      data_min, data_max, data_stp, data_lbl
character(len= 16) :: domain, ch_data
logical :: menue_flag=.false.

interface
  subroutine read_plot_par(n_unit, lbl, val1, val2, val3, val4, val5)
    use sicograph_types
    implicit none
    integer(i4b) :: n_unit
    character(len=256) :: lbl, val1
    character(len=256), optional :: val2, val3, val4, val5
    end subroutine read_plot_par
end interface


!-------- Reading of time-series data --------

call read_core(runname, ndata)

!-------- Menue --------

do

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') '  (1) Temperature-deviation forcing'
   write(6,'(1x,a)') '  (2) Sea-level forcing'
   write(6,'(1x,a)') '  (3) Ice thickness'
   write(6,'(1x,a)') '  (4) Surface velocity'
   write(6,'(1x,a)') '  (5) Basal temperature'
   write(6,'(1x,a)') '  (6) Basal frictional heating'
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' > '
   read (5,*) menue

   if (menue >= 3) then

#if (defined(ANT))
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)') '  (1) Vostok'
      write(6,'(1x,a)') '  (2) Dome A'
      write(6,'(1x,a)') '  (3) Dome C'
      write(6,'(1x,a)') '  (4) Dome F'
      write(6,'(1x,a)') '  (5) Kohnen'
      write(6,'(1x,a)') '  (6) Byrd'
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)',advance='no') ' > '
      read (5,*) n_core
#elif (defined(EMTP2SGE))
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)') '  (1) P1'
      write(6,'(1x,a)') '  (2) P2'
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)',advance='no') ' > '
      read (5,*) n_core
#elif (defined(GRL))
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)') '  (1) GRIP'
      write(6,'(1x,a)') '  (2) GISP2'
      write(6,'(1x,a)') '  (3) Dye3'
      write(6,'(1x,a)') '  (4) Camp Century'
      write(6,'(1x,a)') '  (5) NGRIP'
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)',advance='no') ' > '
      read (5,*) n_core
#elif (defined(HEINO))
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)') '  (1) P1'
      write(6,'(1x,a)') '  (2) P2'
      write(6,'(1x,a)') '  (3) P3'
      write(6,'(1x,a)') '  (4) P4'
      write(6,'(1x,a)') '  (5) P5'
      write(6,'(1x,a)') '  (6) P6'
      write(6,'(1x,a)') '  (7) P7'
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)',advance='no') ' > '
      read (5,*) n_core
#elif (defined(NMARS))
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)') '  (1) NP'
      write(6,'(1x,a)') '  (2) C1'
      write(6,'(1x,a)') '  (3) C2'
      write(6,'(1x,a)') ' '
      write(6,'(1x,a)',advance='no') ' > '
      read (5,*) n_core
#endif

   end if

   if (menue == 1) then

      data    =  D_Ts
      ch_data = 'D_Ts'
      menue_flag = .true.

   else if (menue == 2) then

      data    =  z_sl
      ch_data = 'z_sl'
      menue_flag = .true.

   else if (menue == 3) then

#if (defined(ANT))
      if (n_core == 1) then
         data    =  H_Vo
         ch_data = 'H_Vo'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  H_DA
         ch_data = 'H_DA'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  H_DC
         ch_data = 'H_DC'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  H_DF
         ch_data = 'H_DF'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  H_Ko
         ch_data = 'H_Ko'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  H_By
         ch_data = 'H_By'
         menue_flag = .true.
      end if
#elif (defined(EMTP2SGE))
      if (n_core == 1) then
         data    =  H_P1
         ch_data = 'H_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  H_P2
         ch_data = 'H_P2'
         menue_flag = .true.
      end if
#elif (defined(GRL))
      if (n_core == 1) then
         data    =  HHGR
         ch_data = 'H_GR'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  H_G2
         ch_data = 'H_G2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  H_D3
         ch_data = 'H_D3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  H_CC
         ch_data = 'H_CC'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  H_NG
         ch_data = 'H_NG'
         menue_flag = .true.
      end if
#elif (defined(HEINO))
      if (n_core == 1) then
         data    =  H_P1
         ch_data = 'H_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  H_P2
         ch_data = 'H_P2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  H_P3
         ch_data = 'H_P3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  H_P4
         ch_data = 'H_P4'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  H_P5
         ch_data = 'H_P5'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  H_P6
         ch_data = 'H_P6'
         menue_flag = .true.
      else if (n_core == 7) then
         data    =  H_P7
         ch_data = 'H_P7'
         menue_flag = .true.
      end if
#elif (defined(NMARS))
      if (n_core == 1) then
         data    =  H_NP
         ch_data = 'H_NP'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  H_C1
         ch_data = 'H_C1'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  H_C2
         ch_data = 'H_C2'
         menue_flag = .true.
      end if
#endif

   else if (menue == 4) then

#if (defined(ANT))
      if (n_core == 1) then
         data    =  vs_Vo
         ch_data = 'vs_Vo'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  vs_DA
         ch_data = 'vs_DA'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  vs_DC
         ch_data = 'vs_DC'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  vs_DF
         ch_data = 'vs_DF'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  vs_Ko
         ch_data = 'vs_Ko'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  vs_By
         ch_data = 'vs_By'
         menue_flag = .true.
      end if
#elif (defined(EMTP2SGE))
      if (n_core == 1) then
         data    =  vs_P1
         ch_data = 'vs_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  vs_P2
         ch_data = 'vs_P2'
         menue_flag = .true.
      end if
#elif (defined(GRL))
      if (n_core == 1) then
         data    =  vs_GR
         ch_data = 'vs_GR'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  vs_G2
         ch_data = 'vs_G2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  vs_D3
         ch_data = 'vs_D3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  vs_CC
         ch_data = 'vs_CC'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  vs_NG
         ch_data = 'vs_NG'
         menue_flag = .true.
      end if
#elif (defined(HEINO))
      if (n_core == 1) then
         data    =  vs_P1
         ch_data = 'vs_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  vs_P2
         ch_data = 'vs_P2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  vs_P3
         ch_data = 'vs_P3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  vs_P4
         ch_data = 'vs_P4'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  vs_P5
         ch_data = 'vs_P5'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  vs_P6
         ch_data = 'vs_P6'
         menue_flag = .true.
      else if (n_core == 7) then
         data    =  vs_P7
         ch_data = 'vs_P7'
         menue_flag = .true.
      end if
#elif (defined(NMARS))
      if (n_core == 1) then
         data    =  vs_NP
         ch_data = 'vs_NP'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  vs_C1
         ch_data = 'vs_C1'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  vs_C2
         ch_data = 'vs_C2'
         menue_flag = .true.
      end if
#endif

   else if (menue == 5) then

#if (defined(ANT))
      if (n_core == 1) then
         data    =  Tb_Vo
         ch_data = 'Tb_Vo'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Tb_DA
         ch_data = 'Tb_DA'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Tb_DC
         ch_data = 'Tb_DC'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  Tb_DF
         ch_data = 'Tb_DF'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  Tb_Ko
         ch_data = 'Tb_Ko'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  Tb_By
         ch_data = 'Tb_By'
         menue_flag = .true.
      end if
#elif (defined(EMTP2SGE))
      if (n_core == 1) then
         data    =  Tb_P1
         ch_data = 'Tb_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Tb_P2
         ch_data = 'Tb_P2'
         menue_flag = .true.
      end if
#elif (defined(GRL))
      if (n_core == 1) then
         data    =  Tb_GR
         ch_data = 'Tb_GR'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Tb_G2
         ch_data = 'Tb_G2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Tb_D3
         ch_data = 'Tb_D3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  Tb_CC
         ch_data = 'Tb_CC'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  Tb_NG
         ch_data = 'Tb_NG'
         menue_flag = .true.
      end if
#elif (defined(HEINO))
      if (n_core == 1) then
         data    =  Tb_P1
         ch_data = 'Tb_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Tb_P2
         ch_data = 'Tb_P2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Tb_P3
         ch_data = 'Tb_P3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  Tb_P4
         ch_data = 'Tb_P4'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  Tb_P5
         ch_data = 'Tb_P5'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  Tb_P6
         ch_data = 'Tb_P6'
         menue_flag = .true.
      else if (n_core == 7) then
         data    =  Tb_P7
         ch_data = 'Tb_P7'
         menue_flag = .true.
      end if
#elif (defined(NMARS))
      if (n_core == 1) then
         data    =  Tb_NP
         ch_data = 'Tb_NP'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Tb_C1
         ch_data = 'Tb_C1'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Tb_C2
         ch_data = 'Tb_C2'
         menue_flag = .true.
      end if
#endif

   else if (menue == 6) then

#if (defined(ANT))
      if (n_core == 1) then
         data    =  Rb_Vo
         ch_data = 'Rb_Vo'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Rb_DA
         ch_data = 'Rb_DA'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Rb_DC
         ch_data = 'Rb_DC'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  Rb_DF
         ch_data = 'Rb_DF'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  Rb_Ko
         ch_data = 'Rb_Ko'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  Rb_By
         ch_data = 'Rb_By'
         menue_flag = .true.
      end if
#elif (defined(EMTP2SGE))
      if (n_core == 1) then
         data    =  Rb_P1
         ch_data = 'Rb_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Rb_P2
         ch_data = 'Rb_P2'
         menue_flag = .true.
      end if
#elif (defined(GRL))
      if (n_core == 1) then
         data    =  Rb_GR
         ch_data = 'Rb_GR'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Rb_G2
         ch_data = 'Rb_G2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Rb_D3
         ch_data = 'Rb_D3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  Rb_CC
         ch_data = 'Rb_CC'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  Rb_NG
         ch_data = 'Rb_NG'
         menue_flag = .true.
      end if
#elif (defined(HEINO))
      if (n_core == 1) then
         data    =  Rb_P1
         ch_data = 'Rb_P1'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Rb_P2
         ch_data = 'Rb_P2'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Rb_P3
         ch_data = 'Rb_P3'
         menue_flag = .true.
      else if (n_core == 4) then
         data    =  Rb_P4
         ch_data = 'Rb_P4'
         menue_flag = .true.
      else if (n_core == 5) then
         data    =  Rb_P5
         ch_data = 'Rb_P5'
         menue_flag = .true.
      else if (n_core == 6) then
         data    =  Rb_P6
         ch_data = 'Rb_P6'
         menue_flag = .true.
      else if (n_core == 7) then
         data    =  Rb_P7
         ch_data = 'Rb_P7'
         menue_flag = .true.
      end if
#elif (defined(NMARS))
      if (n_core == 1) then
         data    =  Rb_NP
         ch_data = 'Rb_NP'
         menue_flag = .true.
      else if (n_core == 2) then
         data    =  Rb_C1
         ch_data = 'Rb_C1'
         menue_flag = .true.
      else if (n_core == 3) then
         data    =  Rb_C2
         ch_data = 'Rb_C2'
         menue_flag = .true.
      end if
#endif

   end if

   if (menue_flag) exit

end do

!-------- Writing of data on temporary file --------

filename_with_path = GMT_SCRIPT_PATH//'/tmp/' &
                     //trim(runname)//'_'//trim(ch_data)//'_tmp.ser'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do n=1, ndata
   write(n_unit,'(2(2x,1pe10.3))') time(n), data(n)
end do

close(n_unit, status='keep')

!-------- Reading of plot parameters --------

para_file          = 'time_series_core_'//trim(domain)//'.dat'
filename_with_path = 'parameter_files/'//trim(para_file)

open(n_unit, iostat=ios, file=trim(filename_with_path), status='old')

if (ios == 0) then

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' available.'
   write(6,'(1x,a)') 'Parameters listed in this file used.'

   call read_plot_par(n_unit, time_lbl, time_min, time_max, time_stp)
   do n=1, menue
      call read_plot_par(n_unit, data_lbl, data_min, data_max, data_stp)
   end do

else

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' not available.'
   write(6,'(1x,a)') 'Parameters will be computed automatically.'

   time_lbl = 'none'   ! dummy
   time_min = 'none'   ! dummy
   time_max = 'none'   ! dummy
   time_stp = 'none'   ! dummy

   data_lbl = 'none'   ! dummy
   data_min = 'none'   ! dummy
   data_max = 'none'   ! dummy
   data_stp = 'none'   ! dummy

end if

close(n_unit, status='keep')

!-------- Plotting of data --------

 call  system(GMT_SCRIPT_SHELL//' '//GMT_SCRIPT_PATH &
         //'/time_series.gmt' &
         //' '//trim(domain) &
         //' '//trim(runname)//'_'//trim(ch_data) &
         //' '//trim(ch_data) &
         //' '//trim(adjustl(time_min)) &
         //' '//trim(adjustl(time_max)) &
         //' '//trim(adjustl(time_stp)) &
         //' '//trim(adjustl(data_min)) &
         //' '//trim(adjustl(data_max)) &
         //' '//trim(adjustl(data_stp)) &
         //' '//trim(adjustl(time_lbl)) &
         //' '//trim(adjustl(data_lbl)))

end subroutine time_series_core

#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Scatter plots (simulation results vs. data) for a given time slice.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine scatter_plot(domain, runname)

use sicograph_types
use sicograph_variables

#if (defined(NAG))
use f90_unix_proc
    ! Required for executing shell commands with "call system('...')"
#endif

implicit none

character(len=256), intent(in) :: runname
character(len= 16), intent(in) :: domain

integer(i4b) :: i, j, m, n
integer(i4b) :: menue, menue_read_data
integer(i4b) :: iexit
integer(i4b) :: ios
integer(i1b) :: i_cb_flag, i_cl_flag
integer(i4b) :: ndata, ndata_mean
integer(i4b), parameter :: n_unit=11
integer(i1b), dimension(0:IMAX,0:JMAX) :: maske_gr_sim
real(sp), dimension(0:IMAX,0:JMAX) :: zs_gr_sim, zb_gr_sim, &
                                      zl_gr_sim, zl0_gr_sim, &
                                      H_gr_sim, &
                                      vh_s_gr_sim
real(sp), dimension((IMAX+1)*(JMAX+1)) :: data_obs, data_sim
real(sp) :: data_diff_mean, data_diff_rms
real(sp), parameter :: eps=1.0e-05_sp
character(len=256) :: para_file
character(len=256) :: filename_with_path
character(len=256) :: data_obs_lbl, data_obs_min, data_obs_max, data_obs_stp, &
                      data_sim_lbl, data_sim_min, data_sim_max, data_sim_stp
character(len= 16) :: ch_data, dx, dx_int, cb_flag, cl_flag, dlambda, &
                      dlambda_int, dphi, dphi_int
character(len=  4) :: ergnum
character          :: ch_dummy
logical :: menue_flag=.false.

interface
  subroutine read_plot_par(n_unit, lbl, val1, val2, val3, val4, val5)
    use sicograph_types
    implicit none
    integer(i4b) :: n_unit
    character(len=256) :: lbl, val1
    character(len=256), optional :: val2, val3, val4, val5
    end subroutine read_plot_par
end interface


!-------- Grid spacing --------

#if (GRID==0 || GRID==1)

if (abs(DX-4.0).lt.eps) then
   dx =   '4.0'; dx_int =  '04'
else if (abs(DX-5.0).lt.eps) then
   dx =   '5.0'; dx_int =  '05'
else if (abs(DX-10.0).lt.eps) then
   dx =  '10.0'; dx_int =  '10'
else if (abs(DX-20.0).lt.eps) then
   dx =  '20.0'; dx_int =  '20'
else if (abs(DX-25.0).lt.eps) then
   dx =  '25.0'; dx_int =  '25'
else if (abs(DX-40.0).lt.eps) then
   dx =  '40.0'; dx_int =  '40'
else if (abs(DX-50.0).lt.eps) then
   dx =  '50.0'; dx_int =  '50'
else if (abs(DX-75.0).lt.eps) then
   dx =  '75.0'; dx_int =  '75'
else if (abs(DX-80.0).lt.eps) then
   dx =  '80.0'; dx_int =  '80'
else if (abs(DX-100.0).lt.eps) then
   dx = '100.0'; dx_int = '100'
else if (abs(DX-150.0).lt.eps) then
   dx = '150.0'; dx_int = '150'
else if (abs(DX-160.0).lt.eps) then
   dx = '160.0'; dx_int = '160'
else if (abs(DX-200.0).lt.eps) then
   dx = '200.0'; dx_int = '200'
else
   stop ' >>> scatter_plot: Grid spacing cannot be determined!'
end if

#elif (GRID==2)

if ( (abs(DLAMBDA-(1.0_dp/3.0_dp).lt.eps)) &
     .and. (abs(DPHI-(1.0_dp/3.0_dp).lt.eps)) ) then
   dlambda = '20.0'; dphi = '20.0';  dlambda_int = '20'; dphi_int = '20'
else if ( (abs(DLAMBDA-(1.0_dp/6.0_dp).lt.eps)) &
          .and. (abs(DPHI-(1.0_dp/6.0_dp).lt.eps)) ) then
   dlambda = '10.0'; dphi = '10.0';  dlambda_int = '10'; dphi_int = '10'
else if ( (abs(DLAMBDA-(1.0_dp/12.0_dp).lt.eps)) &
          .and. (abs(DPHI-(1.0_dp/12.0_dp).lt.eps)) ) then
   dlambda = '5.0'; dphi = '5.0';  dlambda_int = '05'; dphi_int = '05'
else
   stop ' >>> scatter_plot: Grid spacing cannot be determined!'
end if

#endif

!-------- Menue --------

do

   write(6,'(1x,a)') ' '
#if (defined(GRL))
   write(6,'(1x,a)') 'Only for the NEGIS area:'
#endif
   write(6,'(1x,a)') '  (2) Ice thickness'
#if (defined(GRL))
   write(6,'(1x,a)') '  (5) Surface velocity'
#endif

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' > '
   read (5,*) menue

   if (menue == 2) then
      menue_read_data = 12
      menue_flag = .true.
   else if (menue == 5) then
      menue_read_data = 20
      menue_flag      = .true.
   end if

   if (menue_flag) exit

end do

!-------- Reading of simulated and observed data --------

call read_erg_nc(runname, ergnum, menue)

zs_gr_sim    = zs_gr
zb_gr_sim    = zb_gr
zl_gr_sim    = zl_gr
zl0_gr_sim   = zl0_gr
maske_gr_sim = maske_gr
H_gr_sim     = H_gr

do i=0, IMAX
do j=0, JMAX
   if (maske_gr_sim(i,j) == 0) then
      vh_s_gr_sim(i,j) = vh_s_gr(i,j)
   else
      vh_s_gr_sim(i,j) = 0.0
   end if
end do
end do

#if (GRID==0 || GRID==1)
   call read_data(domain, dx_int, menue_read_data)
#elif (GRID==2)
   call read_data(domain, dphi_int, menue_read_data)
#endif

!-------- Creation of data vectors --------

n              = 0
ndata          = 0
ndata_mean     = 0
data_diff_mean = 0.0
data_diff_rms  = 0.0

do i=0, IMAX
do j=0, JMAX

   if (maske_gr(i,j) < 2) then

      n = n+1

      if (menue == 2) then

         data_obs(n)   = H_gr(i,j)
         data_sim(n)   = H_gr_sim(i,j)
         ch_data       = 'H'

         if ((data_obs(n) > eps).and.(data_sim(n) > eps)) then
            ndata_mean     = ndata_mean+1
            data_diff_mean = data_diff_mean + (data_sim(n)-data_obs(n))
            data_diff_rms  = data_diff_rms  + (data_sim(n)-data_obs(n))**2
         end if

      else if (menue == 5) then

         data_obs(n)   = vh_s_gr(i,j)
         data_sim(n)   = vh_s_gr_sim(i,j)
         ch_data       = 'vs'

         if ((data_obs(n) > eps).and.(data_sim(n) > eps)) then
            ndata_mean     = ndata_mean+1
            data_diff_mean = data_diff_mean + (data_sim(n)-data_obs(n))
            data_diff_rms  = data_diff_rms  + (data_sim(n)-data_obs(n))**2
         end if

      end if

   end if

end do
end do

ndata = n

data_diff_mean =      data_diff_mean/real(ndata_mean,sp)
data_diff_rms  = sqrt(data_diff_rms /real(ndata_mean,sp))

write(6,'(1x,a)') ' '
write(6,'(1x,a,es11.3,a)') 'Mean of misfit:', data_diff_mean, ' [data units]'
write(6,'(1x,a,es11.3,a)') 'RMS of misfit :', data_diff_rms , ' [data units]'

!-------- Reading of plot parameters --------

para_file          = 'scatter_plot_'//trim(domain)//'.dat'
filename_with_path = 'parameter_files/'//trim(para_file)

open(n_unit, iostat=ios, file=trim(filename_with_path), status='old')

if (ios == 0) then

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' available.'
   write(6,'(1x,a)') 'Parameters listed in this file used.'

   do n=1, menue
      call read_plot_par(n_unit, data_obs_lbl, data_obs_min, data_obs_max, data_obs_stp)
      call read_plot_par(n_unit, data_sim_lbl, data_sim_min, data_sim_max, data_sim_stp)
   end do

else

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)') 'Plot-parameter file '//trim(para_file)//' not available.'
   write(6,'(1x,a)') 'Parameters will be computed automatically.'

   data_obs_lbl = 'none'   ! dummy
   data_obs_min = 'none'   ! dummy
   data_obs_max = 'none'   ! dummy
   data_obs_stp = 'none'   ! dummy
   data_sim_lbl = 'none'   ! dummy
   data_sim_min = 'none'   ! dummy
   data_sim_max = 'none'   ! dummy
   data_sim_stp = 'none'   ! dummy

end if

close(n_unit, status='keep')

!-------- Writing of data on temporary file --------

filename_with_path = GMT_SCRIPT_PATH//'/tmp/' &
                     //trim(runname)//trim(ergnum)//'_' &
                     //trim(ch_data)//'_scatter'//'_tmp.sca'

open(n_unit, iostat=ios, file=trim(filename_with_path), status='replace')

do n=1, ndata
   write(n_unit,'(2(2x,1pe10.3))') data_obs(n), data_sim(n)
end do

close(n_unit, status='keep')

!-------- Plotting of data --------

call system(GMT_SCRIPT_SHELL//' '//GMT_SCRIPT_PATH &
         //'/scatter.gmt' &
         //' '//trim(domain) &
         //' '//trim(runname)//trim(ergnum)//'_'//trim(ch_data)//'_scatter' &
         //' '//trim(ch_data) &
         //' '//trim(adjustl(data_obs_min)) &
         //' '//trim(adjustl(data_obs_max)) &
         //' '//trim(adjustl(data_obs_stp)) &
         //' '//trim(adjustl(data_sim_min)) &
         //' '//trim(adjustl(data_sim_max)) &
         //' '//trim(adjustl(data_sim_stp)) &
         //' '//trim(adjustl(data_obs_lbl)) &
         //' '//trim(adjustl(data_sim_lbl)) )

end subroutine scatter_plot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading of plot parameters.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_plot_par(n_unit, lbl, val1, val2, val3, val4, val5)

use sicograph_types

implicit none
integer(i4b) :: n_unit
character :: check
character(len=256) :: lbl, val1
character(len=256), optional :: val2, val3, val4, val5
character(len=*), parameter :: form1='(a)', form2='(/)', &
          form3='(a24,1x,a15)', form4='(1x,a15)'

do
   read(n_unit,form1,advance='no') check
   if (check == '%') then   ! Comment line
      read(n_unit,form2,advance='no')
   else   ! Data line
      read(n_unit,form3,advance='no') lbl, val1
      if (present(val2)) read(n_unit,form4,advance='no') val2
      if (present(val3)) read(n_unit,form4,advance='no') val3
      if (present(val4)) read(n_unit,form4,advance='no') val4
      if (present(val5)) read(n_unit,form4,advance='no') val5
      read(n_unit,form2,advance='no')
      exit
   end if
end do

end subroutine read_plot_par

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading and processing of data of time-slice files *.erg (unformatted)
!! or *.nc (NetCDF format).
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_erg_nc(runname, ergnum, menue)

use sicograph_types
use sicograph_variables
#if (NETCDF > 1)
use netcdf
#endif

implicit none
integer(i4b) :: i, j, kc, kt, kr, n
integer(i4b) :: ios
integer(i4b) :: menue
integer(i1b) :: flag_3d_output
real(sp) :: BETA
character(len=256) :: filename, filename_with_path
character(len=256) :: runname
character(len=  4) :: ergnum

integer(i1b), dimension(0:IMAX,0:JMAX) :: maske_erg, maske_old_erg
integer(i1b), dimension(0:IMAX,0:JMAX) :: n_cts_erg
integer(i4b), dimension(0:IMAX,0:JMAX) :: kc_cts_erg
real(dp) :: year2sec_erg=0.0_dp, time_erg=0.0_dp, &
            delta_ts_erg=0.0_dp, glac_index_erg=0.0_dp, z_sl_erg=0.0_dp, &
            V_tot_erg=0.0_dp, V_af_erg=0.0_dp, &
            A_grounded_erg=0.0_dp, A_floating_erg=0.0_dp, &
            xi_erg(0:IMAX), eta_erg(0:JMAX), &
            sigma_level_c_erg(0:KCMAX), sigma_level_t_erg(0:KTMAX), &
            sigma_level_r_erg(0:KRMAX)
real(sp) :: H_R_erg
real(sp), dimension(0:IMAX,0:JMAX) :: lambda_erg, phi_erg, &
            lon_erg, lat_erg, &
            temp_s_erg, accum_erg, as_perp_erg, as_perp_apl_erg, &
            q_geo_erg, &
            zs_erg, zm_erg, zb_erg, zl_erg, zl0_erg, &
            H_cold_erg, H_temp_erg, H_erg, &
            Q_bm_erg, Q_tld_erg, &
            am_perp_erg, &
            qx_erg, qy_erg, &
            dzs_dtau_erg, dzm_dtau_erg, dzb_dtau_erg, dzl_dtau_erg, &
            dH_c_dtau_erg, dH_t_dtau_erg, dH_dtau_erg, &
            vx_b_g_erg, vy_b_g_erg, vz_b_erg, vh_b_erg, &
            vx_s_g_erg, vy_s_g_erg, vz_s_erg, vh_s_erg, &
            vx_m_g_erg, vy_m_g_erg,           vh_m_erg, &
            temp_b_erg, temph_b_erg, &
            tau_b_driving_erg, tau_b_drag_erg, &
            p_b_w_erg, H_w_erg, q_gl_g_erg, q_cf_g_erg
real(sp), dimension(0:IMAX,0:JMAX,0:KCMAX) :: vx_c_erg, vy_c_erg, vz_c_erg, &
                                              temp_c_erg, age_c_erg, &
                                              enth_c_erg, omega_c_erg, &
                                              enh_c_erg
real(sp), dimension(0:IMAX,0:JMAX,0:KTMAX) :: vx_t_erg, vy_t_erg, vz_t_erg, &
                                              omega_t_erg, age_t_erg, &
                                              enth_t_erg, &
                                              enh_t_erg
real(sp), dimension(0:IMAX,0:JMAX,0:KRMAX) :: temp_r_erg

#if (DISC>0)   /* Ice discharge parameterisation */
integer(i1b), dimension(0:IMAX,0:JMAX) :: mask_mar_erg
real(sp),     dimension(0:IMAX,0:JMAX) :: dis_perp_erg, &
                                          cst_dist_erg, cos_grad_tc_erg
#endif

#if (NETCDF > 1)
integer(i4b) :: ncid, ncv
!     ncid:      ID of the output file
!     ncv:       Variable ID
#endif

real(dp), parameter :: pi      = 3.141592653589793_dp
real(dp), parameter :: rad2deg = 180.0_dp/pi

real(sp), parameter :: eps = 1.0e-05_sp

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

BETA = 8.70e-04   ! terrestrial Clausius-Clapeyron gradient, in K/m

#elif ( defined(NMARS) \
        || defined(SMARS) )   /* Martian ice sheet */

BETA = 3.30e-04   ! Martian Clausius-Clapeyron gradient, in K/m

#else

stop ' >>> read_erg_nc: No valid domain specified!'

#endif

!-------- Enter name of time-slice file --------

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Name of run: ', trim(runname)
write(6,'(1x,a)',advance='no') &
        'Number of time-slice file (with leading zeros, 4 digits) > '
read (5,'(a)') ergnum
write(6,'(1x,a)') 'Time-slice file contains:'
write(6,'(1x,a)',advance='no') &
        ' (0) only 2-d arrays, (1) full set of 2-d and 3-d arrays > '
read (5,*) flag_3d_output

if ((flag_3d_output /= 0).and.(flag_3d_output /= 1)) flag_3d_output = 1

!-------- Reading of data from time-slice file --------

read_erg_nc_eof_flag = .false.

#if (NETCDF==1)   /* unformatted time-slice file */

if (flag_3d_output == 0) then
   filename = trim(runname)//'_2d_'//trim(ergnum)//'.erg'
else
   filename = trim(runname)//trim(ergnum)//'.erg'
end if

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Expected file name: ', trim(filename)

filename_with_path = trim(OUTPATH)//'/'//trim(filename)

open(10, iostat=ios, file=trim(filename_with_path), status='old', &
         form='unformatted')
if (ios /= 0) stop ' >>> read_erg_nc: Error when opening the file!'

read(10) ch_attr_title
read(10) ch_attr_institution
read(10) ch_attr_source
read(10) ch_attr_history
read(10) ch_attr_references

read(10) year2sec_erg
read(10) time_erg
read(10) delta_ts_erg
read(10) z_sl_erg

read(10) V_tot_erg
read(10) V_af_erg
read(10) A_grounded_erg
read(10) A_floating_erg

read(10) xi_erg
read(10) eta_erg
read(10) sigma_level_c_erg
read(10) sigma_level_t_erg
read(10) sigma_level_r_erg

read(10) lon_erg
read(10) lat_erg
read(10) lambda_erg
read(10) phi_erg
read(10) temp_s_erg
read(10) accum_erg
read(10) as_perp_erg
read(10) as_perp_apl_erg

#if (DISC>0)   /* Ice discharge parameterisation */

read(10) dis_perp_erg
read(10) cst_dist_erg
read(10) cos_grad_tc_erg
read(10) mask_mar_erg

#endif

read(10) q_geo_erg
read(10) maske_erg
read(10) maske_old_erg
read(10) n_cts_erg
read(10) kc_cts_erg
read(10) zs_erg
read(10) zm_erg
read(10) zb_erg
read(10) zl_erg
read(10) zl0_erg
read(10) H_cold_erg
read(10) H_temp_erg
read(10) H_erg
read(10) H_R_erg
read(10) Q_bm_erg
read(10) Q_tld_erg
read(10) am_perp_erg
read(10) qx_erg
read(10) qy_erg
read(10) dzs_dtau_erg
read(10) dzm_dtau_erg
read(10) dzb_dtau_erg
read(10) dzl_dtau_erg
read(10) dH_c_dtau_erg
read(10) dH_t_dtau_erg
read(10) dH_dtau_erg
read(10) vx_b_g_erg
read(10) vy_b_g_erg
read(10) vz_b_erg
read(10) vh_b_erg
read(10) vx_s_g_erg
read(10) vy_s_g_erg
read(10) vz_s_erg
read(10) vh_s_erg
read(10) vx_m_g_erg
read(10) vy_m_g_erg
read(10) vh_m_erg
read(10) temp_b_erg
read(10) temph_b_erg
read(10) tau_b_driving_erg
read(10) tau_b_drag_erg
read(10) p_b_w_erg
read(10) H_w_erg
read(10) q_gl_g_erg
read(10) q_cf_g_erg

if (flag_3d_output == 1) then
   read(10) vx_c_erg
   read(10) vy_c_erg
   read(10) vz_c_erg
   read(10) vx_t_erg
   read(10) vy_t_erg
   read(10) vz_t_erg
   read(10) temp_c_erg
   read(10) omega_t_erg
   read(10) temp_r_erg
   read(10) enth_c_erg
   read(10) enth_t_erg
   read(10) omega_c_erg
   read(10) enh_c_erg
   read(10) enh_t_erg
   read(10) age_c_erg
   read(10) age_t_erg
else
   write(6,'(/a)') ' >>> read_erg_nc: No 3-d fields available.'
   write(6, '(a)') '                  Some variables will be undefined.'
end if

go to 120

110 continue
write(6,'(/a)') ' >>> read_erg_nc: End-of-file condition while reading *.erg.'
write(6, '(a)') '                  Some variables will be undefined.'
read_erg_nc_eof_flag = .true.

120 continue

close(10, status='keep')

#elif (NETCDF==2)   /* time-slice file in NetCDF format */

if (flag_3d_output == 0) then
   filename = trim(runname)//'_2d_'//trim(ergnum)//'.nc'
else
   filename = trim(runname)//trim(ergnum)//'.nc'
end if

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Expected file name: ', trim(filename)

filename_with_path = trim(OUTPATH)//'/'//trim(filename)

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   write(6,'(/a)') ' >>> read_erg_nc: Error when opening the time-slice file!'
   stop
end if

if ( nf90_inq_varid(ncid, 'year2sec', ncv) == nf90_noerr ) then 
   call check( nf90_get_var(ncid, ncv, year2sec_erg) )
else
   year2sec_erg = 0.0_dp
end if

call check( nf90_inq_varid(ncid, 'time', ncv) )
call check( nf90_get_var(ncid, ncv, time_erg) )

if ( nf90_inq_varid(ncid, 'delta_ts', ncv) == nf90_noerr ) then 
   call check( nf90_get_var(ncid, ncv, delta_ts_erg) )
else if ( nf90_inq_varid(ncid, 'glac_index', ncv) == nf90_noerr ) then 
   call check( nf90_get_var(ncid, ncv, glac_index_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Neither variable delta_ts nor glac_index'
   write(6, '(a)') '                  available in read file *.nc.'
   delta_ts_erg   = 0.0_dp
   glac_index_erg = 0.0_dp
end if

call check( nf90_inq_varid(ncid, 'z_sl', ncv) )
call check( nf90_get_var(ncid, ncv, z_sl_erg) )

call check( nf90_inq_varid(ncid, 'V_tot', ncv) )
call check( nf90_get_var(ncid, ncv, V_tot_erg) )

if ( nf90_inq_varid(ncid, 'V_af', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, V_af_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: V_af'
   write(6, '(a)') '                  not available in read file *.nc.'
   V_af_erg = 0.0_dp
end if

call check( nf90_inq_varid(ncid, 'A_grounded', ncv) )
call check( nf90_get_var(ncid, ncv, A_grounded_erg) )

call check( nf90_inq_varid(ncid, 'A_floating', ncv) )
call check( nf90_get_var(ncid, ncv, A_floating_erg) )

call check( nf90_inq_varid(ncid, 'x', ncv) )
call check( nf90_get_var(ncid, ncv, xi_erg) )

call check( nf90_inq_varid(ncid, 'y', ncv) )
call check( nf90_get_var(ncid, ncv, eta_erg) )

call check( nf90_inq_varid(ncid, 'sigma_level_c', ncv) )
call check( nf90_get_var(ncid, ncv, sigma_level_c_erg) )

call check( nf90_inq_varid(ncid, 'sigma_level_t', ncv) )
call check( nf90_get_var(ncid, ncv, sigma_level_t_erg) )

call check( nf90_inq_varid(ncid, 'sigma_level_r', ncv) )
call check( nf90_get_var(ncid, ncv, sigma_level_r_erg) )

call check( nf90_inq_varid(ncid, 'lon', ncv) )
call check( nf90_get_var(ncid, ncv, lon_erg) )

call check( nf90_inq_varid(ncid, 'lat', ncv) )
call check( nf90_get_var(ncid, ncv, lat_erg) )

call check( nf90_inq_varid(ncid, 'lambda', ncv) )
call check( nf90_get_var(ncid, ncv, lambda_erg) )

call check( nf90_inq_varid(ncid, 'phi', ncv) )
call check( nf90_get_var(ncid, ncv, phi_erg) )

call check( nf90_inq_varid(ncid, 'temp_s', ncv) )
call check( nf90_get_var(ncid, ncv, temp_s_erg) )

if ( nf90_inq_varid(ncid, 'prec', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, accum_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable prec'
   write(6, '(a)') '                  not available in read file *.nc.'
   accum_erg = 0.0_sp
end if

call check( nf90_inq_varid(ncid, 'as_perp', ncv) )
call check( nf90_get_var(ncid, ncv, as_perp_erg) )

if ( nf90_inq_varid(ncid, 'as_perp_apl', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, as_perp_apl_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable as_perp_apl'
   write(6, '(a)') '                  not available in read file *.nc.'
   as_perp_apl_erg = 0.0_sp
end if

#if (DISC>0)   /* Ice discharge parameterisation */

if ( nf90_inq_varid(ncid, 'dis_perp', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, dis_perp_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable dis_perp'
   write(6, '(a)') '                  not available in read file *.nc.'
   dis_perp_erg = 0.0_sp
end if

if ( nf90_inq_varid(ncid, 'cst_dist', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, cst_dist_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable cst_dist'
   write(6, '(a)') '                  not available in read file *.nc.'
   cst_dist_erg = 0.0_sp
end if

if ( nf90_inq_varid(ncid, 'cos_grad_tc', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, cos_grad_tc_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable cos_grad_tc'
   write(6, '(a)') '                  not available in read file *.nc.'
   cos_grad_tc_erg = 0.0_sp
end if

if ( nf90_inq_varid(ncid, 'mask_mar', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, mask_mar_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable mask_mar'
   write(6, '(a)') '                  not available in read file *.nc.'
   mask_mar_erg = 0
end if

#endif

if ( nf90_inq_varid(ncid, 'q_geo', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, q_geo_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable q_geo'
   write(6, '(a)') '                  not available in read file *.nc.'
   q_geo_erg = 0.0_sp
end if

call check( nf90_inq_varid(ncid, 'maske', ncv) )
call check( nf90_get_var(ncid, ncv, maske_erg) )

if ( nf90_inq_varid(ncid, 'maske_old', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, maske_old_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable maske_old'
   write(6, '(a)') '                  not available in read file *.nc.'
   maske_old_erg = 0
end if

call check( nf90_inq_varid(ncid, 'n_cts', ncv) )
call check( nf90_get_var(ncid, ncv, n_cts_erg) )

call check( nf90_inq_varid(ncid, 'kc_cts', ncv) )
call check( nf90_get_var(ncid, ncv, kc_cts_erg) )

call check( nf90_inq_varid(ncid, 'zs', ncv) )
call check( nf90_get_var(ncid, ncv, zs_erg) )

call check( nf90_inq_varid(ncid, 'zm', ncv) )
call check( nf90_get_var(ncid, ncv, zm_erg) )

call check( nf90_inq_varid(ncid, 'zb', ncv) )
call check( nf90_get_var(ncid, ncv, zb_erg) )

call check( nf90_inq_varid(ncid, 'zl', ncv) )
call check( nf90_get_var(ncid, ncv, zl_erg) )

if ( nf90_inq_varid(ncid, 'zl0', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, zl0_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable zl0'
   write(6, '(a)') '                  not available in read file *.nc.'
   zl0_erg = 0.0_sp
end if

call check( nf90_inq_varid(ncid, 'H_cold', ncv) )
call check( nf90_get_var(ncid, ncv, H_cold_erg) )

call check( nf90_inq_varid(ncid, 'H_temp', ncv) )
call check( nf90_get_var(ncid, ncv, H_temp_erg) )

call check( nf90_inq_varid(ncid, 'H', ncv) )
call check( nf90_get_var(ncid, ncv, H_erg) )

call check( nf90_inq_varid(ncid, 'H_R', ncv) )
call check( nf90_get_var(ncid, ncv, H_R_erg) )

call check( nf90_inq_varid(ncid, 'Q_bm', ncv) )
call check( nf90_get_var(ncid, ncv, Q_bm_erg) )

call check( nf90_inq_varid(ncid, 'Q_tld', ncv) )
call check( nf90_get_var(ncid, ncv, Q_tld_erg) )

call check( nf90_inq_varid(ncid, 'am_perp', ncv) )
call check( nf90_get_var(ncid, ncv, am_perp_erg) )

call check( nf90_inq_varid(ncid, 'qx', ncv) )
call check( nf90_get_var(ncid, ncv, qx_erg) )

call check( nf90_inq_varid(ncid, 'qy', ncv) )
call check( nf90_get_var(ncid, ncv, qy_erg) )

call check( nf90_inq_varid(ncid, 'dzs_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dzs_dtau_erg) )

call check( nf90_inq_varid(ncid, 'dzm_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dzm_dtau_erg) )

call check( nf90_inq_varid(ncid, 'dzb_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dzb_dtau_erg) )

call check( nf90_inq_varid(ncid, 'dzl_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dzl_dtau_erg) )

call check( nf90_inq_varid(ncid, 'dH_c_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dH_c_dtau_erg) )

call check( nf90_inq_varid(ncid, 'dH_t_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dH_t_dtau_erg) )

call check( nf90_inq_varid(ncid, 'dH_dt', ncv) )
call check( nf90_get_var(ncid, ncv, dH_dtau_erg) )

call check( nf90_inq_varid(ncid, 'vx_b_g', ncv) )
call check( nf90_get_var(ncid, ncv, vx_b_g_erg) )

call check( nf90_inq_varid(ncid, 'vy_b_g', ncv) )
call check( nf90_get_var(ncid, ncv, vy_b_g_erg) )

call check( nf90_inq_varid(ncid, 'vz_b', ncv) )
call check( nf90_get_var(ncid, ncv, vz_b_erg) )

call check( nf90_inq_varid(ncid, 'vh_b', ncv) )
call check( nf90_get_var(ncid, ncv, vh_b_erg) )

call check( nf90_inq_varid(ncid, 'vx_s_g', ncv) )
call check( nf90_get_var(ncid, ncv, vx_s_g_erg) )

call check( nf90_inq_varid(ncid, 'vy_s_g', ncv) )
call check( nf90_get_var(ncid, ncv, vy_s_g_erg) )

call check( nf90_inq_varid(ncid, 'vz_s', ncv) )
call check( nf90_get_var(ncid, ncv, vz_s_erg) )

call check( nf90_inq_varid(ncid, 'vh_s', ncv) )
call check( nf90_get_var(ncid, ncv, vh_s_erg) )

if ( nf90_inq_varid(ncid, 'vx_m_g', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, vx_m_g_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable vx_m_g'
   write(6, '(a)') '                  not available in read file *.nc.'
   vx_m_g_erg = 0.0_sp
end if

if ( nf90_inq_varid(ncid, 'vy_m_g', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, vy_m_g_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable vy_m_g'
   write(6, '(a)') '                  not available in read file *.nc.'
   vy_m_g_erg = 0.0_sp
end if

if ( nf90_inq_varid(ncid, 'vh_m', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, vh_m_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable vh_m'
   write(6, '(a)') '                  not available in read file *.nc.'
   vh_m_erg = 0.0_sp
end if

call check( nf90_inq_varid(ncid, 'temp_b', ncv) )
call check( nf90_get_var(ncid, ncv, temp_b_erg) )

call check( nf90_inq_varid(ncid, 'temph_b', ncv) )
call check( nf90_get_var(ncid, ncv, temph_b_erg) )

if ( nf90_inq_varid(ncid, 'tau_b_driving', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, tau_b_driving_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable tau_b_driving'
   write(6, '(a)') '                  not available in read file *.nc.'
   tau_b_driving_erg = 0.0_sp
end if

if ( nf90_inq_varid(ncid, 'tau_b_drag', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, tau_b_drag_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable tau_b_drag'
   write(6, '(a)') '                  not available in read file *.nc.'
   tau_b_drag_erg = 0.0_sp
end if

call check( nf90_inq_varid(ncid, 'p_b_w', ncv) )
call check( nf90_get_var(ncid, ncv, p_b_w_erg) )

call check( nf90_inq_varid(ncid, 'H_w', ncv) )
call check( nf90_get_var(ncid, ncv, H_w_erg) )

call check( nf90_inq_varid(ncid, 'q_gl_g', ncv) )
call check( nf90_get_var(ncid, ncv, q_gl_g_erg) )

if ( nf90_inq_varid(ncid, 'q_cf_g', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, q_cf_g_erg) )
else
   write(6,'(/a)') ' >>> read_erg_nc: Variable q_cf_g'
   write(6, '(a)') '                  not available in read file *.nc.'
   q_cf_g_erg = 0.0_sp
end if

if (flag_3d_output == 1) then

   call check( nf90_inq_varid(ncid, 'vx_c', ncv) )
   call check( nf90_get_var(ncid, ncv, vx_c_erg) )

   call check( nf90_inq_varid(ncid, 'vy_c', ncv) )
   call check( nf90_get_var(ncid, ncv, vy_c_erg) )

   call check( nf90_inq_varid(ncid, 'vz_c', ncv) )
   call check( nf90_get_var(ncid, ncv, vz_c_erg) )

   call check( nf90_inq_varid(ncid, 'vx_t', ncv) )
   call check( nf90_get_var(ncid, ncv, vx_t_erg) )

   call check( nf90_inq_varid(ncid, 'vy_t', ncv) )
   call check( nf90_get_var(ncid, ncv, vy_t_erg) )

   call check( nf90_inq_varid(ncid, 'vz_t', ncv) )
   call check( nf90_get_var(ncid, ncv, vz_t_erg) )

   call check( nf90_inq_varid(ncid, 'temp_c', ncv) )
   call check( nf90_get_var(ncid, ncv, temp_c_erg) )

   call check( nf90_inq_varid(ncid, 'omega_t', ncv) )
   call check( nf90_get_var(ncid, ncv, omega_t_erg) )

   call check( nf90_inq_varid(ncid, 'temp_r', ncv) )
   call check( nf90_get_var(ncid, ncv, temp_r_erg) )

   call check( nf90_inq_varid(ncid, 'enth_c', ncv) )
   call check( nf90_get_var(ncid, ncv, enth_c_erg) )

   call check( nf90_inq_varid(ncid, 'enth_t', ncv) )
   call check( nf90_get_var(ncid, ncv, enth_t_erg) )

   call check( nf90_inq_varid(ncid, 'omega_c', ncv) )
   call check( nf90_get_var(ncid, ncv, omega_c_erg) )

   if ( nf90_inq_varid(ncid, 'enh_c', ncv) == nf90_noerr ) then
      call check( nf90_get_var(ncid, ncv, enh_c_erg) )
   else
      write(6,'(/a)') ' >>> read_erg_nc: Variable enh_c'
      write(6, '(a)') '                  not available in read file *.nc.'
      enh_c_erg = 0.0_sp
   end if

   if ( nf90_inq_varid(ncid, 'enh_t', ncv) == nf90_noerr ) then
      call check( nf90_get_var(ncid, ncv, enh_t_erg) )
   else
      write(6,'(/a)') ' >>> read_erg_nc: Variable enh_t'
      write(6, '(a)') '                  not available in read file *.nc.'
      enh_t_erg = 0.0_sp
   end if

   call check( nf90_inq_varid(ncid, 'age_c', ncv) )
   call check( nf90_get_var(ncid, ncv, age_c_erg) )

   call check( nf90_inq_varid(ncid, 'age_t', ncv) )
   call check( nf90_get_var(ncid, ncv, age_t_erg) )

else

   write(6,'(/a)') ' >>> read_erg_nc: No 3-d fields available.'
   write(6, '(a)') '                  Some variables will be undefined.'

end if

call check( nf90_close(ncid) )

#endif

!-------- General abbreviations --------

if (DEFORM >= eps) then

   ea = exp(DEFORM) 

   do kc=0, KCMAX
      zeta_c(kc) = real(kc,sp)/real(KCMAX,sp)
      eaz_c(kc)  = exp(DEFORM*zeta_c(kc)) 
      eaz_c_quotient(kc) = (eaz_c(kc)-1.0)/(ea-1.0)
   end do

else

   ea = 1.0

   do kc=0, KCMAX
      zeta_c(kc) = real(kc,sp)/real(KCMAX,sp)
      eaz_c(kc)  = 1.0
      eaz_c_quotient(kc) = zeta_c(kc)
   end do

end if

do kt=0, KTMAX
   zeta_t(kt) = real(kt,sp)/real(KTMAX,sp)
end do

do kr=0, KRMAX
   zeta_r(kr) = real(kr,sp)/real(KRMAX,sp)
end do

!-------- Conversion of read data --------

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

time_gr = real((time_erg*1.0e-03_dp),sp)   ! a -> ka
write(6,'(/4x,a,es13.6,a)') 'time = ', time_erg, ' a'

#else   /* Martian ice sheet */

time_gr = real((time_erg*1.0e-06_dp),sp)   ! a -> Ma
write(6,'(/4x,a,es13.6,a)') 'time = ', time_erg, ' a'

#endif

delta_ts_gr   = real(delta_ts_erg,sp)         ! deg C
glac_index_gr = real(glac_index_erg,sp)       ! 1
z_sl_gr       = real((z_sl_erg*0.001_dp),sp)  ! m -> km

V_tot_gr      = real((V_tot_erg*1.0e-09_dp),sp)       ! m3 -> km3
A_grounded_gr = real((A_grounded_erg*1.0e-06_dp),sp)  ! m2 -> km2
A_floating_gr = real((A_floating_erg*1.0e-06_dp),sp)  ! m2 -> km2

write(6,'(4x,a,es13.6,a)') 'V_tot      = ', V_tot_gr, ' km^3'
write(6,'(4x,a,es13.6,a)') 'A_grounded = ', A_grounded_gr, ' km^2'
write(6,'(4x,a,es13.6,a)') 'A_floating = ', A_floating_gr, ' km^2'

do i=0, IMAX

#if (GRID==0 || GRID==1)
   xi_gr(i) = real((xi_erg(i)*0.001_dp),sp)  ! m -> km
#elif (GRID==2)
   xi_gr(i) = real((xi_erg(i)*rad2deg),sp)   ! rad -> deg
#endif

end do

do j=0, JMAX

#if (GRID==0 || GRID==1)
   eta_gr(j) = real((eta_erg(j)*0.001_dp),sp)  ! m -> km
#elif (GRID==2)
   eta_gr(j) = real((eta_erg(j)*rad2deg),sp)   ! rad -> deg
#endif

end do

H_R_gr = H_R_erg*0.001    ! m -> km

do i=0, IMAX
do j=0, JMAX

   maske_gr(i,j)  = maske_erg(i,j)
   n_cts_gr(i,j)  = n_cts_erg(i,j)
   kc_cts_gr(i,j) = kc_cts_erg(i,j)

   temp_s_gr(i,j)  = temp_s_erg(i,j)            ! deg C
   as_perp_gr(i,j) = as_perp_erg(i,j) *1000.0   ! m/a -> mm/a

   zs_gr(i,j)     = zs_erg(i,j)     *0.001   ! m -> km
   zm_gr(i,j)     = zm_erg(i,j)     *0.001   ! m -> km
   zb_gr(i,j)     = zb_erg(i,j)     *0.001   ! m -> km
   zl_gr(i,j)     = zl_erg(i,j)     *0.001   ! m -> km
   H_cold_gr(i,j) = H_cold_erg(i,j) *0.001   ! m -> km
   H_temp_gr(i,j) = H_temp_erg(i,j) *0.001   ! m -> km
   H_gr(i,j)      = H_erg(i,j)      *0.001   ! m -> km

   temp_b_gr(i,j)  = temp_b_erg(i,j)    ! deg C
   temph_b_gr(i,j) = temph_b_erg(i,j)   ! deg C

   p_b_w_gr(i,j)  = p_b_w_erg(i,j)*0.001   ! Pa -> kPa
   H_w_gr(i,j)    = H_w_erg(i,j)           ! m

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   qx_gr(i,j) = qx_erg(i,j) *0.001  ! m2/a -> km*m/a
   qy_gr(i,j) = qy_erg(i,j) *0.001  ! m2/a -> km*m/a

   q_gl_g_gr(i,j) = q_gl_g_erg(i,j) *0.001  ! m2/a -> km*m/a

   Q_bm_gr(i,j)    = Q_bm_erg(i,j)     ! m/a
   Q_tld_gr(i,j)   = Q_tld_erg(i,j)    ! m/a
   am_perp_gr(i,j) = am_perp_erg(i,j)  ! m/a

   dzs_dtau_gr(i,j)  = dzs_dtau_erg(i,j)   ! m/a
   dzm_dtau_gr(i,j)  = dzm_dtau_erg(i,j)   ! m/a
   dzb_dtau_gr(i,j)  = dzb_dtau_erg(i,j)   ! m/a
   dzl_dtau_gr(i,j)  = dzl_dtau_erg(i,j)   ! m/a
   dH_c_dtau_gr(i,j) = dH_c_dtau_erg(i,j)  ! m/a
   dH_t_dtau_gr(i,j) = dH_t_dtau_erg(i,j)  ! m/a
   dH_dtau_gr(i,j)   = dH_dtau_erg(i,j)    ! m/a

   vx_b_g_gr(i,j) = vx_b_g_erg(i,j)   ! m/a
   vy_b_g_gr(i,j) = vy_b_g_erg(i,j)   ! m/a
   vz_b_gr(i,j)   = vz_b_erg(i,j)     ! m/a
   vh_b_gr(i,j)   = vh_b_erg(i,j)     ! m/a
   vx_s_g_gr(i,j) = vx_s_g_erg(i,j)   ! m/a
   vy_s_g_gr(i,j) = vy_s_g_erg(i,j)   ! m/a
   vz_s_gr(i,j)   = vz_s_erg(i,j)     ! m/a
   vh_s_gr(i,j)   = vh_s_erg(i,j)     ! m/a

#else   /* Martian ice sheet */

   qx_gr(i,j) = qx_erg(i,j)  ! m2/a
   qy_gr(i,j) = qy_erg(i,j)  ! m2/a

   q_gl_g_gr(i,j) = q_gl_g_erg(i,j)  ! m2/a

   Q_bm_gr(i,j)    = Q_bm_erg(i,j)    *1000.0   ! m/a -> mm/a
   Q_tld_gr(i,j)   = Q_tld_erg(i,j)   *1000.0   ! m/a -> mm/a
   am_perp_gr(i,j) = am_perp_erg(i,j) *1000.0   ! m/a -> mm/a

   dzs_dtau_gr(i,j)  = dzs_dtau_erg(i,j)  *1000.0   ! m/a -> mm/a
   dzm_dtau_gr(i,j)  = dzm_dtau_erg(i,j)  *1000.0   ! m/a -> mm/a
   dzb_dtau_gr(i,j)  = dzb_dtau_erg(i,j)  *1000.0   ! m/a -> mm/a
   dzl_dtau_gr(i,j)  = dzl_dtau_erg(i,j)  *1000.0   ! m/a -> mm/a
   dH_c_dtau_gr(i,j) = dH_c_dtau_erg(i,j) *1000.0   ! m/a -> mm/a
   dH_t_dtau_gr(i,j) = dH_t_dtau_erg(i,j) *1000.0   ! m/a -> mm/a
   dH_dtau_gr(i,j)   = dH_dtau_erg(i,j)   *1000.0   ! m/a -> mm/a

   vx_b_g_gr(i,j) = vx_b_g_erg(i,j) *1000.0   ! m/a -> mm/a
   vy_b_g_gr(i,j) = vy_b_g_erg(i,j) *1000.0   ! m/a -> mm/a
   vz_b_gr(i,j)   = vz_b_erg(i,j)   *1000.0   ! m/a -> mm/a
   vh_b_gr(i,j)   = vh_b_erg(i,j)   *1000.0   ! m/a -> mm/a
   vx_s_g_gr(i,j) = vx_s_g_erg(i,j) *1000.0   ! m/a -> mm/a
   vy_s_g_gr(i,j) = vy_s_g_erg(i,j) *1000.0   ! m/a -> mm/a
   vz_s_gr(i,j)   = vz_s_erg(i,j)   *1000.0   ! m/a -> mm/a
   vh_s_gr(i,j)   = vh_s_erg(i,j)   *1000.0   ! m/a -> mm/a

#endif

   do kc=0, KCMAX

#if (CALCMOD==1)
      z_c_gr(i,j,kc) = zm_gr(i,j) + H_cold_gr(i,j)*eaz_c_quotient(kc)   ! km
#elif (CALCMOD==0 || CALCMOD==2)
      z_c_gr(i,j,kc) = zb_gr(i,j) + H_gr(i,j)     *eaz_c_quotient(kc)   ! km
#endif

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

      vx_c_gr(i,j,kc) = vx_c_erg(i,j,kc)   ! m/a
      vy_c_gr(i,j,kc) = vy_c_erg(i,j,kc)   ! m/a
      vz_c_gr(i,j,kc) = vz_c_erg(i,j,kc)   ! m/a

      age_c_gr(i,j,kc) = age_c_erg(i,j,kc) *1.0e-03   ! a -> ka

#else   /* Martian ice sheet */

      vx_c_gr(i,j,kc) = vx_c_erg(i,j,kc) *1000.0   ! m/a -> mm/a
      vy_c_gr(i,j,kc) = vy_c_erg(i,j,kc) *1000.0   ! m/a -> mm/a
      vz_c_gr(i,j,kc) = vz_c_erg(i,j,kc) *1000.0   ! m/a -> mm/a

      age_c_gr(i,j,kc) = age_c_erg(i,j,kc) *1.0e-06   ! a -> Ma

#endif

      temp_c_gr(i,j,kc) = temp_c_erg(i,j,kc)

      enth_c_gr(i,j,kc)  = enth_c_erg(i,j,kc)
      omega_c_gr(i,j,kc) = omega_c_erg(i,j,kc)

      enh_c_gr(i,j,kc)   = enh_c_erg(i,j,kc)

   end do

   do kt=0, KTMAX

#if (CALCMOD==1)
      z_t_gr(i,j,kt) = zb_gr(i,j) + H_temp_gr(i,j)*zeta_t(kt)  ! km
#elif (CALCMOD==0 || CALCMOD==2)
      z_t_gr(i,j,kt) = zb_gr(i,j)                              ! km
#endif

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

      vx_t_gr(i,j,kt) = vx_t_erg(i,j,kt)   ! m/a
      vy_t_gr(i,j,kt) = vy_t_erg(i,j,kt)   ! m/a
      vz_t_gr(i,j,kt) = vz_t_erg(i,j,kt)   ! m/a

      age_t_gr(i,j,kt) = age_t_erg(i,j,kt) *1.0e-03   ! a -> ka

#else   /* Martian ice sheet */

      vx_t_gr(i,j,kt) = vx_t_erg(i,j,kt) *1000.0   ! m/a -> mm/a
      vy_t_gr(i,j,kt) = vy_t_erg(i,j,kt) *1000.0   ! m/a -> mm/a
      vz_t_gr(i,j,kt) = vz_t_erg(i,j,kt) *1000.0   ! m/a -> mm/a

      age_t_gr(i,j,kt) = age_t_erg(i,j,kt) *1.0e-06   ! a -> Ma

#endif

      enth_t_gr(i,j,kt)  = enth_t_erg(i,j,kt)
      omega_t_gr(i,j,kt) = omega_t_erg(i,j,kt)

      enh_t_gr(i,j,kt)   = enh_t_erg(i,j,kt)

   end do

   do kr=0, KRMAX

      z_r_gr(i,j,kr) = zb_gr(i,j) + H_R_gr*(real(kr,sp)/real(KRMAX,sp)-1.0) ! km

      temp_r_gr(i,j,kr) = temp_r_erg(i,j,kr)

   end do

end do
end do

!-------- Computation of temperature relative to pressure melting --------

do i=0, IMAX
do j=0, JMAX
do kc=0, KCMAX

#if (CALCMOD==1)
   temph_c_gr(i,j,kc) = temp_c_gr(i,j,kc) &
      - ( -1000.0*BETA*H_cold_gr(i,j)*(1.0-eaz_c_quotient(kc)) )
          ! Factor 1000.0 required because H_cold_gr(i,j) is in km instead of m
#elif (CALCMOD==0 || CALCMOD==2)
   temph_c_gr(i,j,kc) = temp_c_gr(i,j,kc) &
      - ( -1000.0*BETA*H_gr(i,j)*(1.0-eaz_c_quotient(kc)) )
          ! Factor 1000.0 required because H_gr(i,j) is in km instead of m
#endif

end do
end do
end do

!-------- NetCDF error capturing --------

#if (NETCDF > 1)

contains

   subroutine check(status)
      integer(i4b), intent(in) :: status
      if(status /= nf90_noerr) then 
         write(6,'(1x,a)') trim(nf90_strerror(status))
         stop ' >>> read_erg_nc: Stopped due to netcdf error!'
      end if
   end subroutine check  

#endif

end subroutine read_erg_nc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading of data of time-series files *.ser.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_ser(runname, ndata)

use sicograph_types
use sicograph_variables

implicit none
integer(i4b) :: n, ndata
character(len=256) :: runname
character(len=256) :: filename_with_path
character :: ch_dummy

!-------- Name of time-series file --------

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Name of time-series file: ', trim(runname)//'.ser'

!-------- Reading of data from time-series file --------

filename_with_path = OUTPATH//'/'//trim(runname)//'.ser'

open(10, file=trim(filename_with_path), status='old')

read(10,'(a)') ch_dummy
read(10,'(a)') ch_dummy
read(10,'(a)') ch_dummy
read(10,'(a)') ch_dummy
read(10,'(a)') ch_dummy

do n=1, 1048576

   read(10,*,end=20) time(n), D_Ts(n), z_sl(n)
   read(10,*,end=20) V_tot(n), V_grounded(n), V_floating(n), &
                     A_tot(n), A_grounded(n), A_floating(n)
   read(10,*,end=20) V_sle(n), V_temp(n), A_temp(n)
   read(10,*,end=20) H_max(n), H_t_max(n), zs_max(n), vs_max(n), Tbh_max(n)


   read(10,'(a)',end=20) ch_dummy

#if ( defined(ANT) \
      || defined(GRL) \
      || defined(SCAND) \
      || defined(NHEM) \
      || defined(TIBET) \
      || defined(EMTP2SGE) \
      || defined(HEINO) )   /* terrestrial ice sheet */

   time(n)   = time(n)    * 1.0e-03   ! a -> ka
   V_tot(n)  = V_tot(n)   * 1.0e-15   ! m^3 -> 10^6 km^3
   A_tot(n)  = A_tot(n)   * 1.0e-12   ! m^2 -> 10^6 km^2
   V_temp(n) = V_temp(n)  * 1.0e-12   ! m^3 -> 10^3 km^3
   A_temp(n) = A_temp(n)  * 1.0e-12   ! m^2 -> 10^6 km^2
   H_max(n)  = H_max(n)   * 1.0e-03   ! m -> km
   zs_max(n) = zs_max(n)  * 1.0e-03   ! m -> km

#elif (defined(ASF))   /* Austfonna */

   time(n)   = time(n)    * 1.0e-03   ! a -> ka
   V_tot(n)  = V_tot(n)   * 1.0e-12   ! m^3 -> 10^3 km^3
   A_tot(n)  = A_tot(n)   * 1.0e-09   ! m^2 -> 10^3 km^2
   V_temp(n) = V_temp(n)  * 1.0e-12   ! m^3 -> 10^3 km^3
   A_temp(n) = A_temp(n)  * 1.0e-09   ! m^2 -> 10^3 km^2
   V_sle(n)  = V_sle(n)   * 1.0e+03   ! m -> mm

#else   /* Martian ice sheet */

   time(n)   = time(n)    * 1.0e-06   ! a -> Ma
   V_tot(n)  = V_tot(n)   * 1.0e-15   ! m^3 -> 10^6 km^3
   A_tot(n)  = A_tot(n)   * 1.0e-12   ! m^2 -> 10^6 km^2
   V_temp(n) = V_temp(n)  * 1.0e-12   ! m^3 -> 10^3 km^3
   A_temp(n) = A_temp(n)  * 1.0e-12   ! m^2 -> 10^6 km^2
   H_max(n)  = H_max(n)   * 1.0e-03   ! m -> km
   zs_max(n) = zs_max(n)  * 1.0e-03   ! m -> km
   vs_max(n) = vs_max(n)  * 1.0e+03   ! m/a -> mm/a

#endif

end do

20 close(10, status='keep')

ndata = n-1

end subroutine read_ser

#if (defined(HEINO))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading of data of time-series files *.sed
!! (for the sediment region of the HEINO set-up).
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_sed(runname, ndata)

use sicograph_types
use sicograph_variables

implicit none
integer(i4b) :: n, ndata
character(len=256) :: runname
character(len=256) :: filename_with_path
character :: ch_dummy

!-------- Name of time-series file --------

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Name of time-series file: ', trim(runname)//'.sed'

!-------- Reading of data from time-series file --------

filename_with_path = OUTPATH//'/'//trim(runname)//'.sed'

open(10, file=trim(filename_with_path), status='old')

read(10,'(a)') ch_dummy
read(10,'(a)') ch_dummy
read(10,'(a)') ch_dummy

do n=1, 1048576

   read(10,*,end=20) time(n), D_Ts(n), z_sl(n)
   read(10,*,end=20) H_ave_sed(n), Tbh_ave_sed(n), Atb_sed(n)
   read(10,'(a)',end=20) ch_dummy

   time(n)   = time(n)         * 1.0e-03   ! a -> ka
   H_ave_sed(n) = H_ave_sed(n) * 1.0e-03   ! m -> km
   Atb_sed(n)   = Atb_sed(n)   * 1.0e-12   ! m^2 -> 10^6 km^2

end do

20 close(10, status='keep')

ndata = n-1

end subroutine read_sed

#endif

#if ( defined(ANT) \
      || defined(EMTP2SGE) \
      || defined(GRL) \
      || defined(HEINO) \
      || defined(NMARS) )

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading of data of time-series files *.core (for ice-core locations
!! or selected points).
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_core(runname, ndata)

use sicograph_types
use sicograph_variables

implicit none
integer(i4b) :: n, ndata
character(len=256) :: runname
character(len=256) :: filename_with_path
character :: ch_dummy

!-------- Name of time-series file --------

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Name of time-series file: ', trim(runname)//'.core'

!-------- Reading of data from time-series file --------

filename_with_path = OUTPATH//'/'//trim(runname)//'.core'

open(10, file=trim(filename_with_path), status='old')

do n=1, 6; read(10,'(a)') ch_dummy; end do

do n=1, 1048576

   read(10,*,end=20) time(n), D_Ts(n), z_sl(n)

#if (defined(ANT))

   read(10,*,end=20)  H_Vo(n),  H_DA(n),  H_DC(n),  H_DF(n),  H_Ko(n),  H_By(n)
   read(10,*,end=20) vs_Vo(n), vs_DA(n), vs_DC(n), vs_DF(n), vs_Ko(n), vs_By(n)
   read(10,*,end=20) Tb_Vo(n), Tb_DA(n), Tb_DC(n), Tb_DF(n), Tb_Ko(n), Tb_By(n)
   read(10,*,end=20) Rb_Vo(n), Rb_DA(n), Rb_DC(n), Rb_DF(n), Rb_Ko(n), Rb_By(n)

#elif (defined(EMTP2SGE))

   read(10,*,end=20)  H_P1(n),  H_P2(n)
   read(10,*,end=20) vs_P1(n), vs_P2(n)
   read(10,*,end=20) Tb_P1(n), Tb_P2(n)
   read(10,*,end=20) Rb_P1(n), Rb_P2(n)

#elif (defined(GRL))

   read(10,*,end=20)  HHGR(n),  H_G2(n),  H_D3(n),  H_CC(n),  H_NG(n)
   read(10,*,end=20) vs_GR(n), vs_G2(n), vs_D3(n), vs_CC(n), vs_NG(n)
   read(10,*,end=20) Tb_GR(n), Tb_G2(n), Tb_D3(n), Tb_CC(n), Tb_NG(n)
   read(10,*,end=20) Rb_GR(n), Rb_G2(n), Rb_D3(n), Rb_CC(n), Rb_NG(n)

#elif (defined(HEINO))

   read(10,*,end=20)  H_P1(n),  H_P2(n),  H_P3(n),  H_P4(n),  H_P5(n), &
                      H_P6(n),  H_P7(n)
   read(10,*,end=20) vs_P1(n), vs_P2(n), vs_P3(n), vs_P4(n), vs_P5(n), &
                     vs_P6(n), vs_P7(n)
   read(10,*,end=20) Tb_P1(n), Tb_P2(n), Tb_P3(n), Tb_P4(n), Tb_P5(n), &
                     Tb_P6(n), Tb_P7(n)
   read(10,*,end=20) Rb_P1(n), Rb_P2(n), Rb_P3(n), Rb_P4(n), Rb_P5(n), &
                     Rb_P6(n), Rb_P7(n)

#elif (defined(NMARS))

   read(10,*,end=20)  H_NP(n),  H_C1(n),  H_C2(n)
   read(10,*,end=20) vs_NP(n), vs_C1(n), vs_C2(n)
   read(10,*,end=20) Tb_NP(n), Tb_C1(n), Tb_C2(n)
   read(10,*,end=20) Rb_NP(n), Rb_C1(n), Rb_C2(n)

#endif

   read(10,'(a)',end=20) ch_dummy

end do

20 close(10, status='keep')

ndata = n-1

#if (defined(ANT))

   time = time  * 1.0e-03   ! a -> ka
   H_Vo = H_Vo  * 1.0e-03   ! m -> km
   H_DA = H_DA  * 1.0e-03   ! m -> km
   H_DC = H_DC  * 1.0e-03   ! m -> km
   H_DF = H_DF  * 1.0e-03   ! m -> km
   H_Ko = H_Ko  * 1.0e-03   ! m -> km
   H_By = H_By  * 1.0e-03   ! m -> km

#elif (defined(EMTP2SGE))

   time = time  * 1.0e-03   ! a -> ka
   H_P1 = H_P1  * 1.0e-03   ! m -> km
   H_P2 = H_P2  * 1.0e-03   ! m -> km

#elif (defined(GRL))

   time = time  * 1.0e-03   ! a -> ka
   HHGR = HHGR  * 1.0e-03   ! m -> km
   H_G2 = H_G2  * 1.0e-03   ! m -> km
   H_D3 = H_D3  * 1.0e-03   ! m -> km
   H_CC = H_CC  * 1.0e-03   ! m -> km
   H_NG = H_NG  * 1.0e-03   ! m -> km

#elif (defined(HEINO))

   time = time  * 1.0e-03   ! a -> ka
   H_P1 = H_P1  * 1.0e-03   ! m -> km
   H_P2 = H_P2  * 1.0e-03   ! m -> km
   H_P3 = H_P3  * 1.0e-03   ! m -> km
   H_P4 = H_P4  * 1.0e-03   ! m -> km
   H_P5 = H_P5  * 1.0e-03   ! m -> km
   H_P6 = H_P6  * 1.0e-03   ! m -> km
   H_P7 = H_P7  * 1.0e-03   ! m -> km

#elif (defined(NMARS))

   time  = time  * 1.0e-06   ! a -> Ma
   H_NP  = H_NP  * 1.0e-03   ! m -> km
   H_C1  = H_C1  * 1.0e-03   ! m -> km
   H_C2  = H_C2  * 1.0e-03   ! m -> km
   vs_NP = vs_NP * 1.0e+03   ! m/a -> mm/a
   vs_C1 = vs_C1 * 1.0e+03   ! m/a -> mm/a
   vs_C2 = vs_C2 * 1.0e+03   ! m/a -> mm/a

#endif

end subroutine read_core

#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Reading of the data for the plan-view data plots (present surface and
!! bedrock topography, precipitation rate etc.).
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_data(domain, dx_int, menue)

use sicograph_types
use sicograph_variables

implicit none
integer(i4b) :: i, j, m, n
integer(i4b) :: ios
integer(i4b) :: menue
integer(i4b) :: n_month
real(sp) :: xi0_gr, eta0_gr
real(sp), parameter :: eps=1.0e-05_sp
character(len=256) :: inpath
character(len=256) :: filename_with_path
character(len= 16) :: domain, dx_int
character :: ch_dummy

!-------- Reading of topography data --------

#if ( defined(ANT) \
      || defined(ASF) \
      || defined(GRL) \
      || defined(NHEM) \
      || defined(NMARS) \
      || defined(SCAND) \
      || defined(SMARS) \
      || defined(TIBET) )

filename_with_path = INPATH//'/'//trim(domain)//'/'//ZS_PRESENT_FILE
open(21, iostat=ios, file=trim(filename_with_path), status='old')
if (ios.ne.0) stop ' >>> read_data: Error when opening the zs file!'

#if (defined(ANT))
filename_with_path = INPATH//'/'//trim(domain)//'/'//ZB_PRESENT_FILE
open(25, iostat=ios, file=trim(filename_with_path), status='old')
if (ios.ne.0) stop ' >>> read_data: Error when opening the zb file!'
#endif

filename_with_path = INPATH//'/'//trim(domain)//'/'//ZL_PRESENT_FILE
open(22, iostat=ios, file=trim(filename_with_path), status='old')
if (ios.ne.0) stop ' >>> read_data: Error when opening the zl file!'

filename_with_path = INPATH//'/'//trim(domain)//'/'//ZL0_FILE
open(23, iostat=ios, file=trim(filename_with_path), status='old')
if (ios.ne.0) stop ' >>> read_data: Error when opening the zl0 file!'

filename_with_path = INPATH//'/'//trim(domain)//'/'//MASK_PRESENT_FILE
open(24, iostat=ios, file=trim(filename_with_path), status='old')
if (ios.ne.0) stop ' >>> read_data: Error when opening the mask file!'

do m=1, 6; read(21,'(a)') ch_dummy; end do
#if (defined(ANT))
do m=1, 6; read(25,'(a)') ch_dummy; end do
#endif
do m=1, 6; read(22,'(a)') ch_dummy; end do
do m=1, 6; read(23,'(a)') ch_dummy; end do
do m=1, 6; read(24,'(a)') ch_dummy; end do

do j=JMAX, 0, -1
   read(21,*) (zs_gr(i,j), i=0,IMAX)
#if (defined(ANT))
   read(25,*) (zb_gr(i,j), i=0,IMAX)
#endif
   read(22,*) (zl_gr(i,j), i=0,IMAX)
   read(23,*) (zl0_gr(i,j), i=0,IMAX)
   read(24,2300) (maske_gr(i,j), i=0,IMAX)
end do

close(21, status='keep')
#if (defined(ANT))
close(25, status='keep')
#endif
close(22, status='keep')
close(23, status='keep')
close(24, status='keep')

2300 format(IMAX(i1),i1)

#if (!defined(ANT))
do i=0, IMAX
do j=0, JMAX
   if (maske_gr(i,j) <= 1 ) then
      zb_gr(i,j) = zl_gr(i,j)
   else if (maske_gr(i,j) == 2 ) then
      zs_gr(i,j) = 0.0   ! present-day
      zb_gr(i,j) = 0.0   ! sea level
   else   ! (maske_gr(i,j)==3)
      stop ' >>> read_data: maske_gr(i,j)==3 not allowed for initial topography!'
   end if
end do
end do
#endif

#else

stop ' >>> read_data: Reading of topography data not yet implemented!'

#endif

#if (defined(ANT))

   zs_gr  = zs_gr  *1.0e-03   ! m -> km
   zb_gr  = zb_gr  *1.0e-03   ! m -> km
   zl_gr  = zl_gr  *1.0e-03   ! m -> km
   zl0_gr = zl0_gr *1.0e-03   ! m -> km

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

   z_sl_gr = 0.0   ! present-day sea level

#elif (defined(ASF))

   zs_gr  = zs_gr  *1.0e-03   ! m -> km
   zb_gr  = zb_gr  *1.0e-03   ! m -> km
   zl_gr  = zl_gr  *1.0e-03   ! m -> km
   zl0_gr = zl0_gr *1.0e-03   ! m -> km

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

   z_sl_gr = 0.0   ! present-day sea level

#elif (defined(EMTP2SGE))

   stop ' >>> read_data: Reading of topography data not yet implemented!'

#elif (defined(GRL))

   zs_gr  = zs_gr  *1.0e-03   ! m -> km
   zb_gr  = zb_gr  *1.0e-03   ! m -> km
   zl_gr  = zl_gr  *1.0e-03   ! m -> km
   zl0_gr = zl0_gr *1.0e-03   ! m -> km

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

   z_sl_gr = 0.0   ! present-day sea level

#elif (defined(HEINO))

   stop ' >>> read_data: Reading of topography data not yet implemented!'

#elif (defined(NHEM))

   zs_gr  = zs_gr  *1.0e-03   ! m -> km
   zb_gr  = zb_gr  *1.0e-03   ! m -> km
   zl_gr  = zl_gr  *1.0e-03   ! m -> km
   zl0_gr = zl0_gr *1.0e-03   ! m -> km

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

   z_sl_gr = 0.0   ! present-day sea level

#elif (defined(TIBET))

   zs_gr  = zs_gr  *1.0e-03   ! m -> km
   zb_gr  = zb_gr  *1.0e-03   ! m -> km
   zl_gr  = zl_gr  *1.0e-03   ! m -> km
   zl0_gr = zl0_gr *1.0e-03   ! m -> km

   xi0_gr  = LAMBDA_0         ! Coordinates of the point
   eta0_gr = PHI_0            ! (i,j) = (0,0), in deg

   z_sl_gr = 0.0   ! present-day sea level

#elif (defined(NMARS))

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

#elif (defined(SCAND))

   zs_gr  = zs_gr  *1.0e-03   ! m -> km
   zb_gr  = zb_gr  *1.0e-03   ! m -> km
   zl_gr  = zl_gr  *1.0e-03   ! m -> km
   zl0_gr = zl0_gr *1.0e-03   ! m -> km

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

   z_sl_gr = 0.0   ! present-day sea level

#elif (defined(SMARS))

   xi0_gr  = X0         ! Coordinates of the point
   eta0_gr = Y0         ! (i,j) = (0,0), in km

#endif

!-------- Reading of precipitation data --------

precip_gr = 0.0   ! dummy values

#if ( defined(ANT) \
      || defined(GRL) \
      || defined(NHEM) \
      || defined(SCAND) \
      || defined(TIBET) )

if (menue == 15) then

#if ( defined(ANT) || defined(GRL) )
   filename_with_path = INPATH//'/'//trim(domain)//'/'//PRECIP_PRESENT_FILE
#elif ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' Month: (1) Jan, (2) Feb, ... , (12) Dec > '
   read (5,*) n_month
   filename_with_path = INPATH//'/'//trim(domain)//'/'//PRECIP_MM_PRESENT_FILE
#endif

   open(21, iostat=ios, file=trim(filename_with_path), status='old')
   if (ios /= 0) stop ' >>> read_data: Error when opening the precip file!'

   do m=1, 6; read(21,'(a)') ch_dummy; end do

#if ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   do n=1, n_month-1
      do m=1, 3; read(21,'(a)') ch_dummy; end do
      do j=JMAX, 0, -1
         read(21,*) (precip_gr(i,j), i=0,IMAX)   ! dummies (wrong month)
      end do
   end do
   do m=1, 3; read(21,'(a)') ch_dummy; end do
#endif

   do j=JMAX, 0, -1
      read(21,*) (precip_gr(i,j), i=0,IMAX)
   end do

   close(21, status='keep')

end if

#endif

!-------- Reading of data for LGM anomaly of precipitation rate --------

precip_lgm_anom_gr = 0.0   ! dummy values

#if ( defined(ANT) \
      || defined(GRL) \
      || defined(NHEM) \
      || defined(SCAND) \
      || defined(TIBET) )

if (menue == 16) then

#if ( defined(ANT) || defined(GRL) )
   filename_with_path = INPATH//'/'//trim(domain)//'/'//PRECIP_ANOM_FILE
#elif ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' Month: (1) Jan, (2) Feb, ... , (12) Dec > '
   read (5,*) n_month
   filename_with_path = INPATH//'/'//trim(domain)//'/'//PRECIP_MM_ANOM_FILE
#endif

   open(21, iostat=ios, file=trim(filename_with_path), status='old')
   if (ios /= 0) stop ' >>> read_data: Error when opening the precip anom file!'

   do m=1, 6; read(21,'(a)') ch_dummy; end do

#if ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   do n=1, n_month-1
      do m=1, 3; read(21,'(a)') ch_dummy; end do
      do j=JMAX, 0, -1
         read(21,*) (precip_lgm_anom_gr(i,j), i=0,IMAX)  ! dummies (wrong month)
      end do
   end do
   do m=1, 3; read(21,'(a)') ch_dummy; end do
#endif

   do j=JMAX, 0, -1
      read(21,*) (precip_lgm_anom_gr(i,j), i=0,IMAX)
   end do

   close(21, status='keep')

#if ( defined(ANT) || defined(GRL) )
   precip_lgm_anom_gr = precip_lgm_anom_gr * PRECIP_ANOM_FACT
#elif ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   precip_lgm_anom_gr = precip_lgm_anom_gr * PRECIP_MM_ANOM_FACT
#endif

end if

#endif

!-------- Reading of data for the surface temperature --------

temp_s_gr = 0.0   ! dummy values

#if ( defined(ANT) \
      || defined(GRL) \
      || defined(NHEM) \
      || defined(SCAND) \
      || defined(TIBET) )

if (menue == 17) then

#if (defined(ANT))

   inpath = '/uchi/greve/Documents/sico_data/ant_temp_huybrechts'

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' (1) mean annual, (2) mean January > '
   read (5,*) n_month
   if (n_month == 1) then
      filename_with_path = trim(inpath)//'/'//trim(domain)// &
                           trim(dx_int)//'_temp_ma.dat'
   else if (n_month == 2) then
      filename_with_path = trim(inpath)//'/'//trim(domain)// &
                           trim(dx_int)//'_temp_mj.dat'
   end if

#elif (defined(GRL))

   inpath = '/uchi/greve/Documents/sico_data/grl_temp_ritz'

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' (1) mean annual, (2) mean July > '
   read (5,*) n_month
   if (n_month == 1) then
      filename_with_path = trim(inpath)//'/'//trim(domain)// &
                           trim(dx_int)//'_temp_ma.dat'
   else if (n_month == 2) then
      filename_with_path = trim(inpath)//'/'//trim(domain)// &
                           trim(dx_int)//'_temp_mj.dat'
   end if

#elif ( defined(NHEM) || defined(SCAND) || defined(TIBET) )

   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' Month: (1) Jan, (2) Feb, ... , (12) Dec > '
   read (5,*) n_month
   filename_with_path = INPATH//'/'//trim(domain)//'/'//TEMP_MM_PRESENT_FILE

#endif

   open(21, iostat=ios, file=trim(filename_with_path), status='old')
   if (ios /= 0) stop ' >>> read_data: Error when opening the temp_s file!'

   do m=1, 6; read(21,'(a)') ch_dummy; end do

#if ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   do n=1, n_month-1
      do m=1, 3; read(21,'(a)') ch_dummy; end do
      do j=JMAX, 0, -1
         read(21,*) (temp_s_gr(i,j), i=0,IMAX)   ! dummies (wrong month)
      end do
   end do
   do m=1, 3; read(21,'(a)') ch_dummy; end do
#endif

   do j=JMAX, 0, -1
      read(21,*) (temp_s_gr(i,j), i=0,IMAX)
   end do

   close(21, status='keep')

end if

#endif

!-------- Reading of data for LGM anomaly of surface temperature --------

temp_s_lgm_anom_gr = 0.0   ! dummy values

#if ( defined(ANT) \
      || defined(GRL) \
      || defined(NHEM) \
      || defined(SCAND) \
      || defined(TIBET) )

if (menue == 18) then

#if (defined(ANT))
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' (1) mean annual, (2) mean January > '
   read (5,*) n_month
   if (n_month == 1) then
      filename_with_path = INPATH//'/'//trim(domain)//'/'//TEMP_MA_ANOM_FILE
   else if (n_month == 2) then
      filename_with_path = INPATH//'/'//trim(domain)//'/'//TEMP_MJ_ANOM_FILE
   end if
#elif (defined(GRL))
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' (1) mean annual, (2) mean July > '
   read (5,*) n_month
   if (n_month == 1) then
      filename_with_path = INPATH//'/'//trim(domain)//'/'//TEMP_MA_ANOM_FILE
   else if (n_month == 2) then
      filename_with_path = INPATH//'/'//trim(domain)//'/'//TEMP_MJ_ANOM_FILE
   end if
#elif ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   write(6,'(1x,a)') ' '
   write(6,'(1x,a)',advance='no') ' Month: (1) Jan, (2) Feb, ... , (12) Dec > '
   read (5,*) n_month
   filename_with_path = INPATH//'/'//trim(domain)//'/'//TEMP_MM_ANOM_FILE
#endif

   open(21, iostat=ios, file=trim(filename_with_path), status='old')
   if (ios /= 0) stop ' >>> read_data: Error when opening the temp_s anom file!'

   do m=1, 6; read(21,'(a)') ch_dummy; end do

#if ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   do n=1, n_month-1
      do m=1, 3; read(21,'(a)') ch_dummy; end do
      do j=JMAX, 0, -1
         read(21,*) (temp_s_lgm_anom_gr(i,j), i=0,IMAX)  ! dummies (wrong month)
      end do
   end do
   do m=1, 3; read(21,'(a)') ch_dummy; end do
#endif

   do j=JMAX, 0, -1
      read(21,*) (temp_s_lgm_anom_gr(i,j), i=0,IMAX)
   end do

   close(21, status='keep')

#if ( defined(ANT) || defined(GRL) )
   if (n_month == 1) then
      temp_s_lgm_anom_gr = temp_s_lgm_anom_gr * TEMP_MA_ANOM_FACT
   else if (n_month == 2) then
      temp_s_lgm_anom_gr = temp_s_lgm_anom_gr * TEMP_MJ_ANOM_FACT
   end if
#elif ( defined(NHEM) || defined(SCAND) || defined(TIBET) )
   temp_s_lgm_anom_gr = temp_s_lgm_anom_gr * TEMP_MM_ANOM_FACT
#endif

end if

#endif

!-------- Reading of geothermal-heat-flux data --------

qgeo_gr = 0.0   ! dummy values

#if ( defined(ANT) || defined(GRL) || defined(NMARS) )

if (menue == 19) then

   filename_with_path = INPATH//'/'//trim(domain)//'/'//Q_GEO_FILE

   open(21, iostat=ios, file=trim(filename_with_path), status='old')
   if (ios /= 0) stop ' >>> read_data: Error when opening the qgeo file!'

   do m=1, 6; read(21,'(a)') ch_dummy; end do

   do j=JMAX, 0, -1
      read(21,*) (qgeo_gr(i,j), i=0,IMAX)
   end do

   close(21, status='keep')

end if

#endif

!-------- Reading of surface-velocity data --------

vh_s_gr = 0.0   ! dummy values

#if (defined(GRL))

if (menue == 20) then

   inpath = '/uchi/greve/Documents/sico_data/grl_surfvel_searise_joughin_etal'
   filename_with_path = trim(inpath)//'/'//trim(domain)//'_sr_dev1.2_'// &
                        trim(dx_int)//'_vs.dat'

   open(21, iostat=ios, file=trim(filename_with_path), status='old')
   if (ios /= 0) stop ' >>> read_data: Error when opening the vs file!'

   do m=1, 6; read(21,'(a)') ch_dummy; end do

   do j=JMAX, 0, -1
      read(21,*) (vh_s_gr(i,j), i=0,IMAX)
   end do

   close(21, status='keep')

end if

#endif

!-------- Gridpoints on the x and y axis --------

#if (GRID==0 || GRID==1)

do i=0, IMAX
   xi_gr(i)  = xi0_gr+real(i,sp)*DX
end do

do j=0, JMAX
   eta_gr(j) = eta0_gr+real(j,sp)*DX
end do

#elif (GRID==2)

do i=0, IMAX
   xi_gr(i)  = xi0_gr+real(i,sp)*DLAMBDA
end do

do j=0, JMAX
   eta_gr(j) = eta0_gr+real(j,sp)*DPHI
end do

#endif

!-------- Ice thickness --------

do i=0, IMAX
do j=0, JMAX
   if ((maske_gr(i,j) == 0).or.(maske_gr(i,j) == 3)) then
      H_cold_gr(i,j) = zs_gr(i,j)-zb_gr(i,j)
      H_temp_gr(i,j) = 0.0
   else
      H_cold_gr(i,j) = 0.0
      H_temp_gr(i,j) = 0.0
   end if
   H_gr(i,j) = H_cold_gr(i,j)+H_temp_gr(i,j)
end do
end do

end subroutine read_data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                           End of sicograph.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
