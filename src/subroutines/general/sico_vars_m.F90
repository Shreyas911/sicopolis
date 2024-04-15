!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! Declarations of global variables for SICOPOLIS.
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
!> Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_climatol
   !! Surface-temperature (MAAT) climatology

real(dp), dimension(0:JMAX,0:IMAX) :: smb_climatol
   !! SMB climatology

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_anom
   !! Surface-temperature (MAAT) anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom
   !! SMB anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: dtemp_maat_dz
   !! Surface-temperature (MAAT) vertical gradient

real(dp), dimension(0:JMAX,0:IMAX) :: dsmb_dz
   !! SMB vertical gradient

#endif

real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_in
   !! Prescribed SMB correction

#if (defined(ANT) || defined(GRL)) /* Antarctica or Greenland */

logical :: flag_initmip_asmb
   !! Flag for use of InitMIP SMB anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom_initmip
   !! InitMIP SMB anomaly

#endif

#if (defined(ANT)) /* Antarctica */

logical :: flag_initmip_abmb
   !! Flag for use of InitMIP sub-ice-shelf-melt anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: ab_anom_initmip
   !! InitMIP sub-ice-shelf-melt anomaly

logical :: flag_larmip
   !! Flag for use of LARMIP sub-ice-shelf-melt anomaly

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_larmip_region
   !! LARMIP regions for ice shelf basal melting

real(dp), dimension(0:7) :: ab_anom_larmip
   !! LARMIP sub-ice-shelf-melt anomaly

#endif

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_bm_region
   !! Regions for ice shelf basal melting

#if (FLOATING_ICE_BASAL_MELTING==6)

real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm_present
   !! Equidistant depth points of the
   !! present-day thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm_present
   !! Present-day thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm
   !! Equidistant depth points of the
   !! thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm
   !! Thermal forcing data of the ocean

#endif

#if (defined(ANT) && ICE_SHELF_COLLAPSE_MASK==1) /* Antarctica */

real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
   !! Reference ice thickness for the ice-shelf collapse mask

real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
   !! Ice-shelf collapse mask

#endif

#if (defined(GRL) && RETREAT_MASK==1) /* Greenland */

real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
   !! Reference ice thickness for the retreat mask

real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
   !! Retreat mask

#endif

#if (defined(GRL) && DISC>0) /* Greenland */

integer(i4b) :: disc
integer(i4b) :: n_discharge_call
integer(i4b) :: iter_mar_coa
real(dp)     :: c_dis_0
real(dp)     :: s_dis
real(dp)     :: c_dis_fac
real(dp)     :: T_sub_PD
real(dp)     :: alpha_sub
real(dp)     :: alpha_o
real(dp)     :: m_H
real(dp)     :: m_D
real(dp)     :: r_mar_eff
real(dp)     :: T_sea_freeze
real(dp)     :: dT_glann
real(dp)     :: dT_sub

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_mar
real(dp),     dimension(0:JMAX,0:IMAX) :: c_dis
real(dp),     dimension(0:JMAX,0:IMAX) :: cst_dist
real(dp),     dimension(0:JMAX,0:IMAX) :: cos_grad_tc
real(dp),     dimension(0:JMAX,0:IMAX) :: dis_perp

#if (DISC==2)

integer(i4b) :: glann_time_min
   !! Minimum time of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: glann_time_stp
   !! Time step of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: glann_time_max
   !! Maximum time of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: ndata_glann
   !! Number of data values for the global annual temperature anomaly

integer(i4b) , parameter :: ndata_glann_max = 262143
   !! Maximum allowed value of ndata_glann

real(dp), dimension(0:ndata_glann_max) :: dT_glann_CLIMBER
   !! Data values for the global annual temperature anomaly

#endif

#endif

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */

integer(i4b) :: n_surf
   !! Number of surface points for which time-series data are written 

integer(i4b), parameter :: n_surf_max = 256
   !! Maximum allowed value of n_surf

real(dp), dimension(n_surf_max) :: lambda_surf
   !! Geographical longitude of the prescribed surface points

real(dp), dimension(n_surf_max) :: phi_surf
   !! Geographical latitude of the prescribed surface points

real(dp), dimension(n_surf_max) :: x_surf
   !! Coordinate xi (= x) of the prescribed surface points

real(dp), dimension(n_surf_max) :: y_surf
   !! Coordinate eta (= y) of the prescribed surface points

#endif

end module sico_vars_m
!
