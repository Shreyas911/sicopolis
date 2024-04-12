!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  b o u n d a r y _ m
!
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!!
!!##### Authors
!!
!! Ralf Greve
!!
!!##### License
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
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS. If not, see <https://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!-------------------------------------------------------------------------------
module boundary_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of boundary_m:
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!-------------------------------------------------------------------------------
subroutine boundary(time, dtime, dxi, deta, &
                    delta_ts, glac_index, z_mar)

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myceiling, myfloor
  use sico_maths_m, only: my_erfc
#endif /* Tapenade */

  use netcdf
  use nc_check_m

#if ((MARGIN==2) && (MARINE_ICE_FORMATION==2) && (MARINE_ICE_CALVING==9))
  use calving_m
#endif

#if (defined(GRL) && DISC>0)
  use discharge_workers_m, only: discharge, dis_perp, dT_glann
#endif

  use mask_update_sea_level_m
  use pdd_m

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta

real(dp), intent(out)   :: delta_ts, glac_index, z_mar

! Further return variables
! (defined as global variables in module sico_variables_m):
!
!    accum(j,i), runoff(j,i), as_perp(j,i), calving(j,i), temp_s(j,i)

integer(i4b) :: i, j, n
integer(i4b) :: i_gr, i_kl
integer(i4b) :: n_year_CE
integer(i4b) :: n_hemisphere
integer(i4b) :: ios
real(dp), dimension(0:JMAX,0:IMAX) :: z_sl_old
real(dp) :: z_sl_old_mean
real(dp) :: z_sl_min, t1, t2, t3, t4, t5, t6
real(dp) :: time_in_years
real(dp) :: rho_inv
real(dp) :: time_gr, time_kl
real(dp) :: z_sle_present, z_sle_help

real(dp), dimension(0:JMAX,0:IMAX,0:12) :: precip
real(dp), dimension(0:JMAX,0:IMAX,12)   :: temp_mm
real(dp), dimension(12)                 :: temp_mm_help
real(dp), dimension(0:JMAX,0:IMAX)      :: temp_ma
real(dp), dimension(0:JMAX,0:IMAX)      :: temp_summer
real(dp), dimension(0:JMAX,0:IMAX)      :: temp_ampl
real(dp), dimension(0:JMAX,0:IMAX)      :: temp_diff

real(dp), dimension(0:JMAX,0:IMAX) :: accum_prescribed
real(dp), dimension(0:JMAX,0:IMAX) :: runoff_prescribed

real(dp), dimension(0:JMAX,0:IMAX)    :: gamma_s
real(dp), dimension(0:JMAX,0:IMAX)    :: precip_fact
real(dp), dimension(0:JMAX,0:IMAX,12) :: precip_fact_mm

real(dp), dimension(0:JMAX,0:IMAX) :: s_stat, inv_sqrt2_s_stat, &
                                      beta1, beta2, Pmax, mu
real(dp), dimension(0:JMAX,0:IMAX) :: lambda_lti, temp_lti

real(dp) :: gamma_t
real(dp) :: theta_ma, theta_ma_offset, c_ma, kappa_ma, gamma_ma, &
            theta_ma_1, c_ma_1, gamma_ma_1, &
            theta_ma_2, c_ma_2, gamma_ma_2, &
            theta_ma_3, c_ma_3, gamma_ma_3, &
            zs_sep_1, zs_sep_2, &
            theta_mj, theta_mj_offset, c_mj, kappa_mj, gamma_mj
real(dp) :: sine_factor
real(dp) :: gamma_p, zs_thresh, &
            alpha_p, beta_p, temp_0, alpha_t, beta_t, &
            temp_inv, temp_inv_present, &
            temp_rain, temp_snow, &
            inv_delta_temp_rain_snow, coeff(0:5), frac_solid
real(dp) :: r_aux
character(len=256) :: ch_aux
logical, dimension(0:JMAX,0:IMAX) :: check_point

#if (defined(ANT) && ICE_SHELF_COLLAPSE_MASK==1)
integer(i4b)       :: n_year_CE_isc
character(len= 16) :: ch_year_CE_isc
character(len=256) :: filename_isc
real(dp), dimension(0:IMAX,0:JMAX) :: r_mask_retreat_conv
#endif

#if (defined(GRL) && RETREAT_MASK==1)
integer(i4b)       :: n_year_CE_rtr
character(len= 16) :: ch_year_CE_rtr
character(len=256) :: filename_rtr
real(dp), dimension(0:IMAX,0:JMAX) :: r_mask_retreat_conv
#endif

integer(i4b) :: ncid
!     ncid:      File ID
integer(i4b) :: ncv
!     ncv:       Variable ID
integer(i4b) :: nc3cor(3)
!     nc3cor(3): Corner of a 3-d array
integer(i4b) :: nc3cnt(3)
!     nc3cnt(3): Count of a 3-d array

real(dp), parameter :: &
          inv_twelve = 1.0_dp/12.0_dp, one_third = 1.0_dp/3.0_dp

character(len=64), parameter :: thisroutine = 'boundary'

#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: i_time_in_years
real(dp)     :: temp_val
#endif /* Tapenade */

time_in_years = time*sec2year
n_year_CE     = floor((time_in_years+YEAR_ZERO)+eps_sp_dp)

rho_inv       = 1.0_dp/RHO

!-------- Detect hemisphere --------

if (PHI0 > eps) then
   n_hemisphere =  1   ! northern hemisphere
else if (PHI0 < (-eps)) then
   n_hemisphere = -1   ! southern hemisphere
else
   n_hemisphere =  0
end if

!-------- Initialization of variables --------

z_sl_old      = z_sl
z_sl_old_mean = z_sl_mean

delta_ts   = 0.0_dp
glac_index = 0.0_dp
z_sl       = 0.0_dp
dzsl_dtau  = 0.0_dp
z_mar      = 0.0_dp

!-------- Surface-temperature deviation from present values --------

#if (TSURFACE==1)
delta_ts = DELTA_TS0
!                           ! Steady state with prescribed constant
!                           ! air-temperature deviation
#elif (TSURFACE==3)
delta_ts = SINE_AMPLIT &
           *cos(2.0_dp*pi*time_in_years/SINE_PERIOD) &
           -SINE_AMPLIT
!                           ! Sinusoidal air-temperature forcing
#elif (TSURFACE==4)

!  ------ delta_ts from ice-core record

if (time_in_years < real(grip_time_min,dp)) then
   delta_ts = griptemp(0)
else if (time_in_years < real(grip_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_kl = floor((time_in_years &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
#else /* Tapenade */
   call myfloor((time_in_years &
                -real(grip_time_min,dp))/real(grip_time_stp,dp),i_kl)
#endif /* Normal vs. Tapenade */
   i_kl = max(i_kl, 0)

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_gr = ceiling((time_in_years &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
#else /* Tapenade */
   call myceiling(((time_in_years &
          -real(grip_time_min,dp))/real(grip_time_stp,dp)),i_gr)
#endif /* Normal vs. Tapenade */
   i_gr = min(i_gr, ndata_grip)

   if (i_kl == i_gr) then

      delta_ts = griptemp(i_kl)

   else

      time_kl = (grip_time_min + i_kl*grip_time_stp) *year2sec
      time_gr = (grip_time_min + i_gr*grip_time_stp) *year2sec

      delta_ts = griptemp(i_kl) &
                +(griptemp(i_gr)-griptemp(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the ice-core data

   end if

else
   delta_ts  = griptemp(ndata_grip)
end if

delta_ts = delta_ts * GRIP_TEMP_FACT
!                ! modification by constant factor

!-------- Glacial index --------

#elif (TSURFACE==5)

if (time_in_years < real(gi_time_min,dp)) then
   glac_index = glacial_index(0)
else if (time_in_years < real(gi_time_max,dp)) then

#if !defined(ALLOW_TAPENADE)/* Normal */
   i_kl = floor((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp))
#else /* Tapenade */
   call myfloor(((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp)), i_kl)
#endif /* Normal vs. Tapenade */
   i_kl = max(i_kl, 0)

#if !defined(ALLOW_TAPENADE)/* Normal */
   i_gr = ceiling((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp))
#else /* Tapenade */
   call myceiling(((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp)), i_gr)
#endif /* Normal vs. Tapenade */
   i_gr = min(i_gr, ndata_gi)

   if (i_kl == i_gr) then

      glac_index = glacial_index(i_kl)

   else

      time_kl = (gi_time_min + i_kl*gi_time_stp) *year2sec
      time_gr = (gi_time_min + i_gr*gi_time_stp) *year2sec

      glac_index = glacial_index(i_kl) &
                +(glacial_index(i_gr)-glacial_index(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the glacial-index data

   end if

else
   glac_index  = glacial_index(ndata_gi)
end if

#endif

!-------- Greenland only:
!         Global annual temperature anomaly
!         (for parameterizing the sub-ocean temperature anomaly
!                             for the ice discharge parameterization) --------

#if (defined(GRL))

#if (DISC==1)

dT_glann = 0.0_dp   ! anomaly set to zero

#elif (DISC==2)

!  ------ dT_glann from CLIMBER simulations

if (time_in_years < real(glann_time_min,dp)) then
   dT_glann = dT_glann_CLIMBER(0)
else if (time_in_years < real(glann_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_kl = floor((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp))
   i_gr = min(i_gr, ndata_glann)
#else /* Tapenade */
   call myfloor(((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp)), i_kl)
   i_kl = max(i_kl, 0)

   call myceiling(((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp)), i_gr)
   i_gr = min(i_gr, ndata_glann)
#endif /* Normal vs. Tapenade */

   if (i_kl == i_gr) then

      dT_glann = dT_glann_CLIMBER(i_kl)

   else

      time_kl = (glann_time_min + i_kl*glann_time_stp) *year2sec
      time_gr = (glann_time_min + i_gr*glann_time_stp) *year2sec

      dT_glann = dT_glann_CLIMBER(i_kl) &
                +(dT_glann_CLIMBER(i_gr)-dT_glann_CLIMBER(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the data

   end if

else
   dT_glann  = dT_glann_CLIMBER(ndata_glann)
end if

#endif

#endif

!-------- Sea level --------

#if (SEA_LEVEL==1)

!  ------ Temporally constant sea level

z_sl = Z_SL0

#elif (SEA_LEVEL==3)

!  ------ Time-dependent sea level from data

if (time_in_years < real(specmap_time_min,dp)) then
   z_sl = specmap_zsl(0)
else if (time_in_years < real(specmap_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_kl = floor((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_gr = min(i_gr, ndata_specmap)
#else /* Tapenade */
   call myfloor(((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp)), i_kl)
   i_kl = max(i_kl, 0)

   call myceiling(((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp)), i_gr)
   i_gr = min(i_gr, ndata_specmap)
#endif /* Normal vs. Tapenade */

   if (i_kl == i_gr) then

      z_sl = specmap_zsl(i_kl)

   else

      time_kl = (specmap_time_min + i_kl*specmap_time_stp) *year2sec
      time_gr = (specmap_time_min + i_gr*specmap_time_stp) *year2sec

      z_sl = specmap_zsl(i_kl) &
                +(specmap_zsl(i_gr)-specmap_zsl(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the sea-level data

   end if

else
   z_sl  = specmap_zsl(ndata_specmap)
end if

#else

errormsg = ' >>> boundary: Parameter SEA_LEVEL must be either 1 or 3!'
call error(errormsg)

#endif

!  ------ Mean sea level

z_sl_mean = sum(z_sl*cell_area)/sum(cell_area)

!  ------ Time derivative of the sea level

if ( z_sl_old_mean > -999999.9_dp ) then
   dzsl_dtau = (z_sl-z_sl_old)/dtime
else   ! only dummy value for z_sl_old available
   dzsl_dtau = 0.0_dp
end if

!  ------ Minimum bedrock elevation for extent of marine ice

#if (MARGIN==2)

#if (MARINE_ICE_CALVING==2 || MARINE_ICE_CALVING==3)
z_mar = Z_MAR
#elif (MARINE_ICE_CALVING==4 || MARINE_ICE_CALVING==5)
z_mar = FACT_Z_MAR*z_sl_mean
#elif (MARINE_ICE_CALVING==6 || MARINE_ICE_CALVING==7)
if (z_sl_mean >= -80.0_dp) then
   z_mar = 2.5_dp*z_sl_mean
else
   z_mar = 10.25_dp*(z_sl_mean+80.0_dp)-200.0_dp
end if
z_mar = FACT_Z_MAR*z_mar
#endif

#endif

!  ------ Update of the mask according to the sea level

!    ---- Check all sea and floating-ice points and their direct
!         neighbours

do i=0, IMAX
do j=0, JMAX
   check_point(j,i) = .false.
end do
end do

do i=1, IMAX-1
do j=1, JMAX-1
   if (mask(j,i) >= 2) then
      check_point(j  ,i  ) = .true.
      check_point(j  ,i+1) = .true.
      check_point(j  ,i-1) = .true.
      check_point(j+1,i  ) = .true.
      check_point(j-1,i  ) = .true.
   end if
end do
end do

do i=1, IMAX-1
do j=1, JMAX-1
   if (check_point(j,i)) then
      mask_new(j,i) = mask_update_sea_level(i, j)
   end if
end do
end do

!    ---- Assign new values of the mask

do i=1, IMAX-1
do j=1, JMAX-1
   if (check_point(j,i)) then
      mask(j,i) = mask_new(j,i)
   end if
end do
end do

!-------- Surface air temperature --------

#if (defined(ANT) || defined(GRL)) /* Antarctica or Greenland */

gamma_t = 0.0_dp ! topographic lapse rate already part of the parameterization
                 ! for the present-day temperature

#else /* other than Antarctica or Greenland */

#if (defined(TOPO_LAPSE_RATE))
gamma_t = TOPO_LAPSE_RATE * 1.0e-03_dp
             ! topographic lapse rate; K/km -> K/m
#else
gamma_t = 6.5e-03_dp
             ! topographic lapse rate (K/m), default value
#endif

#endif

!  ------ Parameterized present-day temperature (only Antarctica, Greenland)

#if (TSURFACE<=5)

#if (defined(ANT) || defined(GRL)) /* Antarctica or Greenland */

#if (defined(ANT)) /* Antarctica */

#if (TEMP_PRESENT_PARA==1)   /* Parameterization by Fortuin and Oerlemans */
                             !  (1990) for the whole ice sheet

#if (defined(TEMP_PRESENT_OFFSET))
theta_ma_offset = TEMP_PRESENT_OFFSET
theta_mj_offset = TEMP_PRESENT_OFFSET
#else
theta_ma_offset = 0.0_dp
theta_mj_offset = 0.0_dp
#endif

theta_ma = 34.461_dp + theta_ma_offset
gamma_ma = -9.14e-03_dp
c_ma     = -0.688_dp

theta_mj = 16.81_dp + theta_mj_offset
gamma_mj = -6.92e-03_dp
c_mj     = -0.27973_dp

#elif (TEMP_PRESENT_PARA==2)   /* Parameterization by Fortuin and Oerlemans */
                               !  (1990), separately for three different
                               !  elevation ranges

#if (defined(TEMP_PRESENT_OFFSET))
theta_ma_offset = TEMP_PRESENT_OFFSET
theta_mj_offset = TEMP_PRESENT_OFFSET
#else
theta_ma_offset = 0.0_dp
theta_mj_offset = 0.0_dp
#endif

zs_sep_1   =  200.0_dp
zs_sep_2   = 1500.0_dp

theta_ma_1 =  49.642_dp + theta_ma_offset
gamma_ma_1 =   0.0_dp
c_ma_1     =  -0.943_dp

theta_ma_2 =  36.689_dp + theta_ma_offset
gamma_ma_2 =  -5.102e-03_dp
c_ma_2     =  -0.725_dp

theta_ma_3 =   7.405_dp + theta_ma_offset
gamma_ma_3 = -14.285e-03_dp
c_ma_3     =  -0.180_dp

theta_mj   =  16.81_dp + theta_mj_offset
gamma_mj   =  -6.92e-03_dp
c_mj       =  -0.27937_dp

#else

errormsg = ' >>> boundary: Parameter TEMP_PRESENT_PARA must be either 1 or 2!'
call error(errormsg)

#endif

do i=0, IMAX
do j=0, JMAX

!    ---- Present-day mean-annual temperature

#if (TEMP_PRESENT_PARA==1)
   temp_ma_present(j,i) = theta_ma + gamma_ma*zs(j,i) &
                                   + c_ma*abs(phi(j,i))*rad2deg
#elif (TEMP_PRESENT_PARA==2)
   if ( zs(j,i) <= zs_sep_1 ) then
      temp_ma_present(j,i) = theta_ma_1 + gamma_ma_1*zs(j,i) &
                                        + c_ma_1*abs(phi(j,i))*rad2deg
   else if ( zs(j,i) <= zs_sep_2 ) then
      temp_ma_present(j,i) = theta_ma_2 + gamma_ma_2*zs(j,i) &
                                        + c_ma_2*abs(phi(j,i))*rad2deg
   else
      temp_ma_present(j,i) = theta_ma_3 + gamma_ma_3*zs(j,i) &
                                        + c_ma_3*abs(phi(j,i))*rad2deg
   end if
#endif

!    ---- Present-day mean-January (mid-summer) temperature

   temp_mj_present(j,i) = theta_mj + gamma_mj*zs(j,i) &
                                   + c_mj*abs(phi(j,i))*rad2deg

end do
end do

#elif (defined(GRL)) /* Greenland */

#if (TEMP_PRESENT_PARA==1)   /* Parameterization by Ritz et al. (1997) */

#if (defined(TEMP_PRESENT_OFFSET))
theta_ma_offset = TEMP_PRESENT_OFFSET
theta_mj_offset = TEMP_PRESENT_OFFSET
#else
theta_ma_offset = 0.0_dp
theta_mj_offset = 0.0_dp
#endif

theta_ma = 49.13_dp + theta_ma_offset
gamma_ma = -7.992e-03_dp
c_ma     = -0.7576_dp
kappa_ma =  0.0_dp

theta_mj = 30.38_dp + theta_mj_offset
gamma_mj = -6.277e-03_dp
c_mj     = -0.3262_dp
kappa_mj =  0.0_dp

#elif (TEMP_PRESENT_PARA==2)   /* Parameterization by Fausto et al. (2009) */

#if (defined(TEMP_PRESENT_OFFSET))
theta_ma_offset = TEMP_PRESENT_OFFSET
theta_mj_offset = TEMP_PRESENT_OFFSET
#else
theta_ma_offset = 0.0_dp
theta_mj_offset = 0.0_dp
#endif

theta_ma = 41.83_dp + theta_ma_offset
gamma_ma = -6.309e-03_dp
c_ma     = -0.7189_dp
kappa_ma = -0.0672_dp

theta_mj = 14.70_dp + theta_mj_offset
gamma_mj = -5.426e-03_dp
c_mj     = -0.1585_dp
kappa_mj = -0.0518_dp

#else

errormsg = ' >>> boundary: Parameter TEMP_PRESENT_PARA must be either 1 or 2!'
call error(errormsg)

#endif

do i=0, IMAX
do j=0, JMAX

!  ------ Present-day mean-annual air temperature

   temp_ma_present(j,i) = theta_ma &
                  + gamma_ma*zs(j,i) &
                  + c_ma*phi(j,i)*rad2deg &
                  + kappa_ma*(modulo(lambda(j,i)+pi,2.0_dp*pi)-pi)*rad2deg
                              ! western longitudes counted negatively

!  ------ Present-day mean-July (mid-summer) air temperature

   temp_mj_present(j,i) = theta_mj &
                    + gamma_mj*zs(j,i) &
                    + c_mj*phi(j,i)*rad2deg &
                    + kappa_mj*(modulo(lambda(j,i)+pi,2.0_dp*pi)-pi)*rad2deg
                                ! western longitudes counted negatively

end do
end do

#endif

!    ---- Amplitude of the annual cycle

do i=0, IMAX
do j=0, JMAX

   temp_ampl(j,i) = temp_mj_present(j,i) - temp_ma_present(j,i)

   if (temp_ampl(j,i) < eps) then
      temp_ampl(j,i) = eps   ! Correction of amplitude, if required
   end if

end do
end do

!    ---- Monthly temperatures

do n=1, 12   ! month counter

   if (n_hemisphere == 1) then   ! northern hemisphere
      sine_factor = sin((real(n,dp)-4.0_dp)*pi/6.0_dp)
   else if (n_hemisphere == -1) then   ! southern hemisphere
      sine_factor = -sin((real(n,dp)-4.0_dp)*pi/6.0_dp)
   else
      sine_factor = 0.0_dp
   end if

   do i=0, IMAX
   do j=0, JMAX
      temp_present(j,i,n) = temp_ma_present(j,i) + sine_factor*temp_ampl(j,i)
   end do
   end do

end do

#endif /* Antarctica or Greenland */

#endif /* (TSURFACE<=5) */

do i=0, IMAX
do j=0, JMAX

#if (TSURFACE<=4)

!  ------ Correction of present monthly temperature with elevation changes
!         and temperature deviation delta_ts

   temp_diff(j,i) = gamma_t*(zs_ref(j,i)-zs(j,i)) + delta_ts

   do n=1, 12   ! month counter
      temp_mm(j,i,n) = temp_present(j,i,n) + temp_diff(j,i)
   end do

#elif (TSURFACE==5)

!  ------ Correction of present monthly temperature with LGM anomaly and
!         glacial index as well as elevation changes

   temp_diff(j,i) = gamma_t*(zs_ref(j,i)-zs(j,i))

   do n=1, 12   ! month counter
      temp_mm(j,i,n) = temp_present(j,i,n) &
                       + glac_index*temp_lgm_anom(j,i,n) &
                       + temp_diff(j,i)
   end do

#endif

!  ------ Mean annual air temperature

#if (TSURFACE<=5)

   temp_ma(j,i) = 0.0_dp   ! initialization value

   do n=1, 12   ! month counter
      temp_ma(j,i) = temp_ma(j,i) + temp_mm(j,i,n)*inv_twelve
   end do

#elif (TSURFACE==6)

   ! ...

#endif

end do
end do

!  ------ Save mean-annual air temperature

temp_maat = temp_ma

!-------- Accumulation-ablation function as_perp --------

#if (ACCSURFACE<=3)

#if (ELEV_DESERT==1)
gamma_p   = GAMMA_P*1.0e-03_dp   ! Precipitation lapse rate
                                 ! for elevation desertification, in m^(-1)
zs_thresh = ZS_THRESH            ! Elevation threshold, in m
#endif

#if (ACCSURFACE==2 || ACCSURFACE==3)
gamma_s = GAMMA_S
#else
gamma_s = 0.0_dp
#endif

#elif (ACCSURFACE==4)

alpha_p =  22.47_dp
beta_p  =   0.046_dp
temp_0  = 273.15_dp
alpha_t =   0.67_dp
beta_t  =  88.9_dp

#endif

#if (ACCSURFACE<=5)

#if (SOLID_PRECIP==1)     /* Marsiat (1994) */

temp_rain =    7.0_dp   ! Threshold monthly mean temperature for
                        ! precipitation = 100% rain, in degC
temp_snow =  -10.0_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% snow, in degC

inv_delta_temp_rain_snow = 1.0_dp/(temp_rain-temp_snow)

#elif (SOLID_PRECIP==2)   /* Bales et al. (2009) */

temp_rain =    7.2_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% rain, in degC
temp_snow =  -11.6_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% snow, in degC

coeff(0) =  5.4714e-01_dp   ! Coefficients
coeff(1) = -9.1603e-02_dp   ! of
coeff(2) = -3.314e-03_dp    ! the
coeff(3) =  4.66e-04_dp     ! fifth-order
coeff(4) =  3.8e-05_dp      ! polynomial
coeff(5) =  6.0e-07_dp      ! fit

#elif (SOLID_PRECIP==3)   /* Huybrechts and de Wolde (1999) */

temp_rain = 2.0_dp      ! Threshold instantaneous temperature for &
                        ! precipitation = 100% rain, in degC
temp_snow = temp_rain   ! Threshold instantaneous temperature for &
                        ! precipitation = 100% snow, in degC

#if (defined(S_STAT_0))
s_stat = S_STAT_0    ! Standard deviation of the air termperature
                     ! (same parameter as in the PDD model)
#else
errormsg = ' >>> boundary: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif

do i=0, IMAX
do j=0, JMAX
   inv_sqrt2_s_stat(j,i) = 1.0_dp/(sqrt(2.0_dp)*s_stat(j,i))
end do
end do

#endif

#endif

#if (ABLSURFACE==1 || ABLSURFACE==2)

#if (defined(S_STAT_0) && defined(BETA1_0) && defined(BETA2_0) && defined(PMAX_0) && defined(MU_0))
s_stat = S_STAT_0
beta1  = BETA1_0  *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*degC) -> (m IE)/(s*degC)
beta2  = BETA2_0  *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*degC) -> (m IE)/(s*degC)
Pmax   = PMAX_0
mu     = MU_0     *(1000.0_dp*86400.0_dp)*(RHO/RHO_W)
                           ! (d*degC)/(mm WE) -> (s*degC)/(m IE)
#else
errormsg = ' >>> boundary: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif

#elif (ABLSURFACE==3)

lambda_lti = LAMBDA_LTI  *(0.001_dp*sec2year)*(RHO_W/RHO)
                         ! (mm WE)/(a*degC) -> (m IE)/(s*degC)
temp_lti   = TEMP_LTI
mu         = 0.0_dp      ! no superimposed ice considered

#endif

do i=0, IMAX
do j=0, JMAX

!  ------ Accumulation

#if (ACCSURFACE<=3)

!    ---- Elevation desertification of precipitation

#if (ELEV_DESERT==0)

   precip_fact(j,i) = 1.0_dp   ! no elevation desertification

#elif (ELEV_DESERT==1)

   if (zs_ref(j,i) < zs_thresh) then
      precip_fact(j,i) &
         = exp(gamma_p*(max(zs(j,i),zs_thresh)-zs_thresh))
   else
      precip_fact(j,i) &
         = exp(gamma_p*(max(zs(j,i),zs_thresh)-zs_ref(j,i)))
   end if

#else
   errormsg = ' >>> boundary: Parameter ELEV_DESERT must be either 0 or 1!'
   call error(errormsg)
#endif

   do n=1, 12   ! month counter
      precip(j,i,n) = precip_present(j,i,n)*precip_fact(j,i)
   end do

#endif

!    ---- Precipitation change related to changing climate

#if (ACCSURFACE==1)
   precip_fact(j,i) = ACCFACT
#elif (ACCSURFACE==2)
   precip_fact(j,i) = 1.0_dp + gamma_s(j,i)*delta_ts
#elif (ACCSURFACE==3)
   precip_fact(j,i) = exp(gamma_s(j,i)*delta_ts)
#endif

#if (ACCSURFACE<=3)

   precip(j,i,0) = 0.0_dp   ! initialization value for mean annual precip

   do n=1, 12   ! month counter
      precip(j,i,n) = precip(j,i,n)*precip_fact(j,i)
                                                  ! monthly precip
      precip(j,i,0) = precip(j,i,0) + precip(j,i,n)*inv_twelve
                                                  ! mean annual precip
   end do

#elif (ACCSURFACE==4)

   precip(j,i,0) = 0.0_dp   ! initialization value for mean annual precip

   temp_inv         = alpha_t * (temp_ma(j,i)+temp_0)         + beta_t   ! in K
   temp_inv_present = alpha_t * (temp_ma_present(j,i)+temp_0) + beta_t   ! in K

   precip_fact(j,i) = exp(alpha_p*(temp_0/temp_inv_present-temp_0/temp_inv)) &
                         *(temp_inv_present/temp_inv)**2 &
                         *(1.0_dp+beta_p*(temp_inv-temp_inv_present))

   do n=1, 12   ! month counter
      precip(j,i,n) = precip_present(j,i,n)*precip_fact(j,i)
                                                  ! monthly precip
      precip(j,i,0) = precip(j,i,0) + precip(j,i,n)*inv_twelve
                                                  ! mean annual precip
   end do

#elif (ACCSURFACE==5)

   precip(j,i,0) = 0.0_dp   ! initialization value for mean annual precip

   do n=1, 12   ! month counter

#if (PRECIP_ANOM_INTERPOL==1)
      precip_fact_mm(j,i,n) = (1.0_dp-glac_index) &
                                 + glac_index*precip_lgm_anom(j,i,n)
                              ! interpolation with a linear function
#elif (PRECIP_ANOM_INTERPOL==2)
      precip_fact_mm(j,i,n) = exp(-glac_index*gamma_precip_lgm_anom(j,i,n))
                              ! interpolation with an exponential function
#endif

      precip(j,i,n) = precip_present(j,i,n)*precip_fact_mm(j,i,n)
                                                      ! monthly precip
      precip(j,i,0) = precip(j,i,0) + precip(j,i,n)*inv_twelve
                                                      ! mean annual precip
   end do

#endif

!    ---- Annual accumulation, snowfall and rainfall rates

#if (ACCSURFACE<=5)

   accum(j,i) = precip(j,i,0)

   snowfall(j,i) = 0.0_dp   ! initialization value

   do n=1, 12   ! month counter

#if (SOLID_PRECIP==1)     /* Marsiat (1994) */

      if (temp_mm(j,i,n) >= temp_rain) then
         frac_solid = 0.0_dp
      else if (temp_mm(j,i,n) <= temp_snow) then
         frac_solid = 1.0_dp
      else
         frac_solid = (temp_rain-temp_mm(j,i,n))*inv_delta_temp_rain_snow
      end if

#elif (SOLID_PRECIP==2)   /* Bales et al. (2009) */

      if (temp_mm(j,i,n) >= temp_rain) then
         frac_solid = 0.0_dp
      else if (temp_mm(j,i,n) <= temp_snow) then
         frac_solid = 1.0_dp
      else
         frac_solid = coeff(0) + temp_mm(j,i,n) * ( coeff(1) &
                               + temp_mm(j,i,n) * ( coeff(2) &
                               + temp_mm(j,i,n) * ( coeff(3) &
                               + temp_mm(j,i,n) * ( coeff(4) &
                               + temp_mm(j,i,n) *   coeff(5) ) ) ) )
               ! evaluation of 5th-order polynomial by Horner scheme
      end if

#elif (SOLID_PRECIP==3)   /* Huybrechts and de Wolde (1999) */

#if !defined(ALLOW_TAPENADE) /* Normal */
      frac_solid = 1.0_dp &
                   - 0.5_dp &
                       *erfc((temp_rain-temp_mm(j,i,n))*inv_sqrt2_s_stat(j,i))
#else /* Tapenade */
      call my_erfc((temp_rain-temp_mm(j,i,n))*inv_sqrt2_s_stat(j,i), temp_val)
      frac_solid = 1.0_dp - 0.5_dp*temp_val
#endif /* Normal vs. Tapenade */

#endif

      snowfall(j,i) = snowfall(j,i) + precip(j,i,n)*frac_solid*inv_twelve

   end do

   rainfall(j,i) = precip(j,i,0) - snowfall(j,i)

   if (snowfall(j,i) < 0.0_dp) snowfall(j,i) = 0.0_dp   ! correction of
   if (rainfall(j,i) < 0.0_dp) rainfall(j,i) = 0.0_dp   ! negative values

#endif

!  ------ Ablation

!    ---- Runoff

#if (ABLSURFACE==1 || ABLSURFACE==2)

!      -- Temperature excess ET

   do n=1, 12   ! month counter
      temp_mm_help(n) = temp_mm(j,i,n)
   end do

   call pdd(temp_mm_help, s_stat(j,i), ET(j,i))

!      -- Formation rate of superimposed ice (melt_star), melt rate (melt)
!         and runoff rate (runoff)

#if (ABLSURFACE==1)

   if ((beta1(j,i)*ET(j,i)) <= (Pmax(j,i)*snowfall(j,i))) then
      melt_star(j,i) = beta1(j,i)*ET(j,i)
      melt(j,i)      = 0.0_dp
      runoff(j,i)    = melt(j,i)+rainfall(j,i)
   else
      melt_star(j,i) = Pmax(j,i)*snowfall(j,i)
      melt(j,i)      = beta2(j,i)*(ET(j,i)-melt_star(j,i)/beta1(j,i))
      runoff(j,i)    = melt(j,i)+rainfall(j,i)
   end if

#elif (ABLSURFACE==2)

   if ( rainfall(j,i) <= (Pmax(j,i)*snowfall(j,i)) ) then

      if ( (rainfall(j,i)+beta1(j,i)*ET(j,i)) &
                 <= (Pmax(j,i)*snowfall(j,i)) ) then
         melt_star(j,i) = rainfall(j,i)+beta1(j,i)*ET(j,i)
         melt(j,i)      = 0.0_dp
         runoff(j,i)    = melt(j,i)
      else
         melt_star(j,i) = Pmax(j,i)*snowfall(j,i)
         melt(j,i)      = beta2(j,i) &
                          *(ET(j,i)-(melt_star(j,i)-rainfall(j,i))/beta1(j,i))
         runoff(j,i)    = melt(j,i)
      end if

   else

      melt_star(j,i) = Pmax(j,i)*snowfall(j,i)
      melt(j,i)      = beta2(j,i)*ET(j,i)
      runoff(j,i)    = melt(j,i) + rainfall(j,i)-Pmax(j,i)*snowfall(j,i)

   end if

#endif

#elif (ABLSURFACE==3)

   if (n_hemisphere == 1) then   ! northern hemisphere
      temp_summer(j,i) = one_third &
                            *(temp_mm(j,i,6)+temp_mm(j,i,7)+temp_mm(j,i,8))
   else if (n_hemisphere == -1) then   ! southern hemisphere
      temp_summer(j,i) = one_third &
                            *(temp_mm(j,i,12)+temp_mm(j,i,1)+temp_mm(j,i,2))
   else
      temp_summer(j,i) = temp_ma(j,i)
   end if

   melt_star(j,i) = 0.0_dp   ! no superimposed ice considered
   melt(j,i)      = lambda_lti(j,i) &
                       *max((temp_summer(j,i)-temp_lti(j,i)), 0.0_dp)
   runoff(j,i)    = melt(j,i) + rainfall(j,i)

#endif

end do
end do

!  ------ SMB = precipitation minus runoff

#if (ACCSURFACE<=5 && ABLSURFACE<=5)

as_perp = accum - runoff

#endif

!    ---- Prescribed SMB correction

smb_corr_prescribed = smb_corr_in

#if (defined(ANT) || defined(GRL))

if (flag_initmip_asmb) then   ! Correction for ISMIP6 InitMIP

#if defined(ALLOW_TAPENADE) /* Tapenade */
   call myfloor(time_in_years, i_time_in_years)
#endif /* Tapenade */

   if ((time_in_years > 0.0_dp).and.(time_in_years <= 40.0_dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
      smb_corr_prescribed = smb_corr_prescribed &
                              + 0.025_dp*floor(time_in_years) * smb_anom_initmip
#else /* Tapenade */
      smb_corr_prescribed = smb_corr_prescribed &
                              + 0.025_dp*(i_time_in_years) * smb_anom_initmip
#endif /* Normal vs. Tapenade */

   else if (time_in_years > 40.0_dp) then

      smb_corr_prescribed = smb_corr_prescribed + smb_anom_initmip

   end if

end if

#endif

as_perp = as_perp + smb_corr_prescribed

accum_prescribed  =  max(smb_corr_prescribed, 0.0_dp)
runoff_prescribed = -min(smb_corr_prescribed, 0.0_dp)

accum  = accum  + accum_prescribed
runoff = runoff + runoff_prescribed

!  ------ Ice-surface temperature (10-m firn temperature) temp_s,
!         including empirical firn-warming correction due to
!         refreezing meltwater when superimposed ice is formed

#if (TSURFACE<=5)

do i=0, IMAX
do j=0, JMAX
   if (melt_star(j,i) >= melt(j,i)) then
      temp_s(j,i) = temp_ma(j,i) + mu(j,i)*(melt_star(j,i)-melt(j,i))
   else
      temp_s(j,i) = temp_ma(j,i)
   end if
end do
end do

#endif

do i=0, IMAX
do j=0, JMAX
   if (temp_s(j,i) > -0.001_dp) then
      temp_s(j,i) = -0.001_dp   ! Cut-off of positive air temperature
   end if
end do
end do

!-------- Calving --------

calving = 0.0_dp   ! Initialization

#if ((MARGIN==2) && (MARINE_ICE_FORMATION==2) && (MARINE_ICE_CALVING==9))

call calving_underwater_ice()

#endif

#if (defined(GRL) && DISC>0) /* Ice discharge parameterization for Greenland */

call discharge(dxi, deta)
calving = calving + dis_perp

#endif

!-------- Antarctica only: Ice-shelf collapse mask --------

#if (defined(ANT) && ICE_SHELF_COLLAPSE_MASK==1)

n_year_CE_isc = n_year_CE

if (n_year_CE_isc > ICE_SHELF_COLLAPSE_MASK_TIME_MAX) &
                       n_year_CE_isc = ICE_SHELF_COLLAPSE_MASK_TIME_MAX

if (firstcall%boundary) r_mask_retreat = 1.0_dp   ! initialization

if (n_year_CE_isc < ICE_SHELF_COLLAPSE_MASK_TIME_MIN) then

   r_mask_retreat = 1.0_dp

else if (n_year_CE_isc /= n_year_CE_isc_save) then

   write(ch_year_CE_isc, '(i0)') n_year_CE_isc

   filename_isc = trim(ICE_SHELF_COLLAPSE_MASK_DIR)//'/'// &
                  trim(ICE_SHELF_COLLAPSE_MASK_FILES)// &
                  trim(ch_year_CE_isc)//'.nc'

   ios = nf90_open(trim(filename_isc), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> boundary: Error when opening the file for the' &
               //                end_of_line &
               //'               ice-shelf collapse mask!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'mask', ncv), thisroutine )
   call check( nf90_get_var(ncid, ncv, r_mask_retreat_conv), thisroutine )

   call check( nf90_close(ncid), thisroutine )

   do i=0, IMAX
   do j=0, JMAX
      r_mask_retreat(j,i) = 1.0_dp - r_mask_retreat_conv(i,j)
                                    ! swap 0 <-> 1
      r_mask_retreat(j,i) = max(min(r_mask_retreat(j,i), 1.0_dp), 0.0_dp)
                                    ! constrain to interval [0,1]
   end do
   end do

end if

n_year_CE_isc_save = n_year_CE_isc

#endif

!-------- Greenland only: Retreat mask due to oceanic forcing --------

#if (defined(GRL) && RETREAT_MASK==1)

n_year_CE_rtr = n_year_CE

if (n_year_CE_rtr > RETREAT_MASK_TIME_MAX) &
                       n_year_CE_rtr = RETREAT_MASK_TIME_MAX

if (firstcall%boundary) r_mask_retreat = 1.0_dp   ! initialization

if (n_year_CE_rtr < RETREAT_MASK_TIME_MIN) then

   r_mask_retreat = 1.0_dp

else if (n_year_CE_rtr /= n_year_CE_rtr_save) then

   write(ch_year_CE_rtr, '(i0)') n_year_CE_rtr

   filename_rtr = trim(RETREAT_MASK_DIR)//'/'// &
                  trim(RETREAT_MASK_FILES)// &
                  trim(ch_year_CE_rtr)//'.nc'

   ios = nf90_open(trim(filename_rtr), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> boundary: Error when opening the file for the' &
               //                end_of_line &
               //'               retreat mask due to oceanic forcing!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'sftgif', ncv), thisroutine )
   call check( nf90_get_var(ncid, ncv, r_mask_retreat_conv), thisroutine )

   call check( nf90_close(ncid), thisroutine )

   do i=0, IMAX
   do j=0, JMAX
      r_mask_retreat(j,i) = max(min(r_mask_retreat_conv(i,j), 1.0_dp), 0.0_dp)
                                    ! constrain to interval [0,1]
   end do
   end do

end if

n_year_CE_rtr_save = n_year_CE_rtr

#endif

!-------- First-call flag --------

if (firstcall%boundary) firstcall%boundary = .false.

end subroutine boundary

!-------------------------------------------------------------------------------

end module boundary_m
!
