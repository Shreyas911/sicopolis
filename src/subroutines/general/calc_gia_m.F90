!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ g i a _ m
!
!> @file
!!
!! Computation of the glacial isostatic adjustment of the lithosphere surface.
!!
!! @section Copyright
!!
!! Copyright 2009-2021 Ralf Greve, Sascha Knell
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

!-------------------------------------------------------------------------------
!> Computation of the glacial isostatic adjustment of the lithosphere surface.
!<------------------------------------------------------------------------------

module calc_gia_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  private
  public :: calc_gia

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_gia_m:
!! Computation of the glacial isostatic adjustment of the lithosphere surface.
!<------------------------------------------------------------------------------
subroutine calc_gia(time, z_sl, dtime, dxi, deta, itercount, iter_wss)

implicit none

integer(i4b), intent(in) :: itercount, iter_wss
real(dp),     intent(in) :: time
real(dp),     intent(in) :: z_sl
real(dp),     intent(in) :: dtime, dxi, deta

integer(i4b) :: i, j
real(dp) :: tldt_inv(0:JMAX,0:IMAX)
real(dp) :: dtime_inv
real(dp) :: rho_rhoa_ratio, rhosw_rhoa_ratio
real(dp) :: time_dimless
real(dp) :: target_topo_tau_factor, target_topo_tau

!-------- Term abbreviations --------

do i=0, IMAX
do j=0, JMAX
   tldt_inv(j,i) = 1.0_dp/(time_lag_asth(j,i)+dtime)
end do
end do

rho_rhoa_ratio   = RHO/RHO_A
rhosw_rhoa_ratio = RHO_SW/RHO_A

dtime_inv = 1.0_dp/dtime

!-------- Computation of zl_new and its time derivative --------

#if (REBOUND==0)

do i=0, IMAX
do j=0, JMAX
   zl_new(j,i)   = zl(j,i)
   dzl_dtau(j,i) = (zl_new(j,i)-zl(j,i))*dtime_inv
end do
end do

#elif (REBOUND==1)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1_i1b) then   ! grounded ice or ice-free land

      zl_new(j,i) = tldt_inv(j,i)*( time_lag_asth(j,i)*zl(j,i) &
                    + dtime*(zl0(j,i) &
                             -FRAC_LLRA*rho_rhoa_ratio*(H_c(j,i)+H_t(j,i))) )

   else   ! (mask(j,i) >= 2_i1b)

      zl_new(j,i) = tldt_inv(j,i)*( time_lag_asth(j,i)*zl(j,i) &
                    + dtime*(zl0(j,i) &
                             -FRAC_LLRA*rhosw_rhoa_ratio*z_sl) )

                   ! Water load relative to the present sea-level stand (0 m)
                   ! -> can be positive or negative

   end if

   dzl_dtau(j,i) = (zl_new(j,i)-zl(j,i))*dtime_inv

end do
end do

#elif (REBOUND==2)

if ((mod(itercount, iter_wss)==0).or.(itercount==1)) then
   write(6, fmt='(10x,a)') 'Computation of wss'
   call calc_elra(z_sl, dxi, deta)   ! Update of the steady-state displacement
                                     ! of the lithosphere (wss)
end if

do i=0, IMAX
do j=0, JMAX
   zl_new(j,i) = tldt_inv(j,i) &
                 * ( time_lag_asth(j,i)*zl(j,i) + dtime*(zl0(j,i)-wss(j,i)) )
   dzl_dtau(j,i) = (zl_new(j,i)-zl(j,i))*dtime_inv
end do
end do

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2)

if (time >= time_target_topo_final) then

   zl_new   = zl_target
   dzl_dtau = 0.0_dp

else if (time > time_target_topo_init) then

   time_dimless =  (time                  -time_target_topo_init) &
                  /(time_target_topo_final-time_target_topo_init)

   time_dimless = max(min(time_dimless, 1.0_dp), 0.0_dp)
                                ! constrain to interval [0,1]

   target_topo_tau_factor = time_dimless*time_dimless*time_dimless &
                               *(10.0_dp + time_dimless &
                                           *(-15.0_dp+6.0_dp*time_dimless))
                                ! make transition smooth (quintic function)

   target_topo_tau_factor = 1.0_dp-target_topo_tau_factor

   target_topo_tau = target_topo_tau_0 * target_topo_tau_factor

   zl_new =   ( target_topo_tau*zl_new + dtime*zl_target ) &
            / ( target_topo_tau        + dtime )

   dzl_dtau = (zl_new-zl)*dtime_inv

end if

#elif (THK_EVOL==3)

target_topo_tau = target_topo_tau_0

zl_new =   ( target_topo_tau*zl_new + dtime*zl_target ) &
         / ( target_topo_tau        + dtime )

dzl_dtau = (zl_new-zl)*dtime_inv

#endif

!======== Generation of an isostatically relaxed
!         lithosphere surface topography. Not to be used regularly! ========

#if (defined(EXEC_MAKE_ZL0))

call make_zl0()

errormsg = ' >>> calc_gia: Routine make_zl0 successfully completed,' &
         //         end_of_line &
         //'        isostatically relaxed lithosphere' &
         //         end_of_line &
         //'        topography written on file' &
         //         end_of_line &
         //'        (in directory specified by OUT_PATH).' &
         //         end_of_line &
         //'        Execution of SICOPOLIS stopped.'
call error(errormsg)   ! actually not an error,
                       ! just a regular stop with an info message

#endif

end subroutine calc_gia

!-------------------------------------------------------------------------------
!> Computation of the isostatic steady-state displacement of the lithosphere
!! (wss) for the ELRA model.
!<------------------------------------------------------------------------------
subroutine calc_elra(z_sl, dxi, deta)

#if defined(ALLOW_OPENAD) /* OpenAD */
  use ctrl_m, only: myfloor
#endif /* OpenAD */

implicit none

real(dp), intent(in) :: z_sl, dxi, deta

integer(i4b) :: i, j, ir, jr, il, jl, n
integer(i4b) :: ir_max, jr_max, min_imax_jmax
integer(i4b) :: il_begin, il_end, jl_begin, jl_end
real(dp)                              :: rho_g, rho_sw_g, rho_a_g_inv
real(dp)                              :: dxi_inv, deta_inv
real(dp)                              :: kei_r_incr_inv
real(dp), dimension(0:JMAX,0:IMAX)    :: l_r, l_r_inv, fac_wss
real(dp)                              :: ra_max

#if !defined(ALLOW_OPENAD) /* Normal */
real(dp), allocatable, dimension(:,:) :: f_0
#else /* OpenAD */
real(dp), dimension(-JMAX:2*JMAX,-IMAX:2*IMAX) :: f_0
#endif /* Normal vs. OpenAD */

real(dp), parameter :: r_infl = 8.0_dp
                       ! Radius of non-locality influence of the elastic
                       ! lithosphere plate (in units of l_r)

!-------- Initialisations --------

rho_g       = RHO*G
rho_sw_g    = RHO_SW*G
rho_a_g_inv = 1.0_dp/(RHO_A*G)

dxi_inv   = 1.0_dp/dxi
deta_inv  = 1.0_dp/deta

kei_r_incr_inv = 1.0_dp/kei_r_incr

do i=0, IMAX
do j=0, JMAX
   l_r(j,i)     = (flex_rig_lith(j,i)*rho_a_g_inv)**0.25_dp
   l_r_inv(j,i) = 1.0_dp/l_r(j,i)
   fac_wss(j,i) = l_r(j,i)*l_r(j,i)/(2.0_dp*pi*flex_rig_lith(j,i))
end do
end do

ra_max  = r_infl*maxval(l_r)   ! Radius of non-locality influence (in m)

#if !defined(ALLOW_OPENAD) /* Normal */

ir_max = floor(ra_max*dxi_inv)
jr_max = floor(ra_max*deta_inv)
         ! Radius of non-locality influence (in grid points)

#else /* OpenAD */

call myfloor(ra_max*dxi_inv, ir_max)
call myfloor(ra_max*deta_inv, jr_max)
         ! Radius of non-locality influence (in grid points)

if (ir_max.gt.IMAX) then
  write(*,*) "ERROR: unhandled case ir_max greater than IMAX for f_0 dimensions"
end if
if (jr_max.gt.JMAX) then
  write(*,*) "ERROR: unhandled case jr_max greater than JMAX for f_0 dimensions"
end if

#endif /* Normal vs. OpenAD */

min_imax_jmax = min(IMAX, JMAX)
ir_max        = min(ir_max, min_imax_jmax)
jr_max        = min(jr_max, min_imax_jmax)
                ! ir_max and jr_max constrained by size of domain

il_begin = -ir_max
il_end   = ir_max+IMAX
jl_begin = -jr_max
jl_end   = jr_max+JMAX

!-------- Ice/water load --------

#if !defined(ALLOW_OPENAD) /* Normal */
allocate(f_0(jl_begin:jl_end, il_begin:il_end))
#endif /* Normal */

do il=il_begin, il_end
do jl=jl_begin, jl_end

   i = min(max(il, 0), IMAX)
   j = min(max(jl, 0), JMAX)

   if (mask(j,i)==0_i1b) then
      f_0(jl,il) = rho_g * area(j,i) * (H_c(j,i) + H_t(j,i))
   else if (mask(j,i)==1_i1b) then
      f_0(jl,il) = 0.0_dp
   else   ! (mask(j,i)>=2_i1b)
      f_0(jl,il) = rho_sw_g * area(j,i) * z_sl
                   ! Water load relative to the present sea-level stand (0 m)
                   ! -> can be positive or negative
   end if

end do
end do

!-------- Computation of the steady-state displacement of the lithosphere
!         (wss, positive downward) --------

do i=0, IMAX
do j=0, JMAX

   wss(j,i) = 0.0_dp

   do ir=-ir_max, ir_max
   do jr=-jr_max, jr_max

#if (GRID==0)
      n = nint( dist_dxdy(jr,ir)*l_r_inv(j,i)*kei_r_incr_inv )
#elif (GRID==1)
      n = nint( sq_g11_g(j,i)*dist_dxdy(jr,ir)*l_r_inv(j,i)*kei_r_incr_inv )
              ! The correction with g11 only is OK because g11=g22 for this case
#elif (GRID==2)
      n = nint( dist_dxdy(jr,ir)*l_r_inv(j,i)*kei_r_incr_inv )
              ! The distortion due to g11 and g22 is already accounted for in
              ! the variable dist_dxdy(jr,ir)
#endif

      wss(j,i) = wss(j,i) - fac_wss(j,i)*f_0(jr+j,ir+i)*kei(n)

   end do
   end do

end do
end do

#if !defined(ALLOW_OPENAD) /* Normal */
deallocate (f_0)
#endif /* Normal */

end subroutine calc_elra

!-------------------------------------------------------------------------------
!> Generation of an isostatically relaxed lithosphere surface topography
!! for either the rigid lithosphere, the local lithosphere or the elastic
!! lithosphere model (depending on the setting of the parameter REBOUND).
!! This routine is not to be used regularly, and it is only executed if the
!! parameter EXEC_MAKE_ZL0 is defined in the header file.
!<------------------------------------------------------------------------------
subroutine make_zl0()

#if defined(ALLOW_OPENAD) /* OpenAD */
  use ctrl_m, only: myceiling 
#endif /* OpenAD */

implicit none

integer(i4b)                       :: i, j, m, n, i_f, j_f, n_filter
integer(i4b)                       :: ios
real(dp), dimension(0:JMAX,0:IMAX) :: zl0_raw, zl0_smoothed
real(dp)                           :: dx
real(dp)                           :: rho_ratio
real(dp)                           :: filter_width, sigma_filter
real(dp)                           :: dist, weigh, sum_weigh
character(len=  8)                 :: ch_resolution
character(len= 64)                 :: ch_model
character(len=256)                 :: filename, filename_with_path

!-------- Checking of the grid --------

#if (GRID==0 || GRID==1)

dx = real(DX, dp)   ! horizontal grid spacing (in km)

#else

dx = 0.0_dp   ! dummy value

errormsg = ' >>> make_zl0: Routine works only for GRID 0 or 1!'
call error(errormsg)

#endif

!-------- Computation of the raw (unsmoothed) topography --------

rho_ratio = RHO/RHO_A

#if (REBOUND==0)

zl0_raw = zl   ! rigid lithosphere

#elif (REBOUND==1)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 0_i1b) then
      zl0_raw(j,i) = zl(j,i) + rho_ratio*(H_c(j,i)+H_t(j,i))
   else
      zl0_raw(j,i) = zl(j,i)
   end if            ! local lithosphere
   
end do
end do

#elif (REBOUND==2)

zl0_raw = zl + wss   ! elastic lithosphere

#else

errormsg = ' >>> make_zl0: REBOUND must be either 0, 1 or 2!'
call error(errormsg)

#endif

!-------- Smoothing of the topography (Gaussian filter) --------

#if (defined(ZL0_FILTER_WIDTH))
  filter_width = real(ZL0_FILTER_WIDTH, dp)
                      ! filter width (half span of filtered area), in km
#else
  filter_width = dx   ! default value: one grid spacing, in km
#endif

sigma_filter = filter_width/dx   ! half span of filtered area,
                                 ! in grid points

#if !defined(ALLOW_OPENAD) /* Normal */
n_filter = ceiling(2.0_dp*sigma_filter)
#else /* OpenAD */
call myceiling(2.0_dp*sigma_filter, n_filter)
#endif /* Normal vs. OpenAD */

n_filter = max(n_filter, 5)

zl0_smoothed = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   sum_weigh = 0.0_dp

   do m=-n_filter, n_filter
   do n=-n_filter, n_filter

      i_f = i+m
      j_f = j+n

      if (i_f <    0) i_f =    0
      if (i_f > IMAX) i_f = IMAX

      if (j_f <    0) j_f =    0
      if (j_f > JMAX) j_f = JMAX

      dist      = sqrt(real(m,dp)**2+real(n,dp)**2)
      weigh     = exp(-(dist/sigma_filter)**2)
      sum_weigh = sum_weigh + weigh

      zl0_smoothed(j,i) = zl0_smoothed(j,i) + weigh*zl0_raw(j_f,i_f)

   end do
   end do

   zl0_smoothed(j,i) = zl0_smoothed(j,i)/sum_weigh

end do
end do

!-------- Writing on file --------

if ( abs(dx-nint(dx)) < eps ) then
   write(ch_resolution, fmt='(i8)') nint(dx)
else
   write(ch_resolution, fmt='(f8.2)') dx
end if

#if !defined(ALLOW_OPENAD) /* Normal */
ch_resolution = adjustl(ch_resolution)
#endif /* Normal */

#if (REBOUND==0)
ch_model = 'rigid_lithosphere'
#elif (REBOUND==1)
ch_model = 'local_lithosphere'
#elif (REBOUND==2)
ch_model = 'elastic_lithosphere'
#endif

filename           = trim(ch_domain_short)//'_'//trim(ch_resolution)
filename           = trim(filename)//'_zl0_'//trim(ch_model)//'.dat'
filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

open(23, iostat=ios, file=trim(filename_with_path), recl=rcl1, status='replace')

if (ios /= 0) then
   errormsg = ' >>> make_zl0: Error when opening the new zl0 file!'
   call error(errormsg)
end if

write(23,'(a)') &
'%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(23,'(a)') '% '//trim(ch_domain_long)//':'
write(23,'(a)') &
'% Relaxed lithosphere surface topography without ice load, in m AMSL.'
write(23,'(a)') &
'% Horizontal resolution '//trim(ch_resolution)//' km.'
write(23,'(a,i3,a,i3,a,i3,a,i3,a)') &
'% ', JMAX+1, ' records [j = ', &
JMAX, ' (-1) 0] with ', IMAX+1, &
' values [i = 0 (1) ', IMAX, '] in each record.'
write(23,'(a)') &
'%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

do j=JMAX, 0, -1
   do i=0, IMAX-1
      zl0_smoothed(j,i) = max(zl0_smoothed(j,i), no_value_neg_2)
      write(23,'(f8.1)', advance='no') zl0_smoothed(j,i)
   end do
   i=IMAX
      zl0_smoothed(j,i) = max(zl0_smoothed(j,i), no_value_neg_2)
      write(23,'(f8.1)') zl0_smoothed(j,i)
end do

close(23, status='keep')

end subroutine make_zl0

!-------------------------------------------------------------------------------

end module calc_gia_m
!
