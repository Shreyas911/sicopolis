!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t h k _ m
!
!> @file
!!
!! Computation of the ice thickness.
!!
!! @section Copyright
!!
!! Copyright 2009-2021 Ralf Greve, Reinhard Calov, Tatsuru Sato
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
!> Computation of the ice thickness.
!<------------------------------------------------------------------------------
module calc_thk_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX), save :: H_new_flow
  real(dp), dimension(0:JMAX,0:IMAX), save :: mb_source
  logical,                            save :: flag_solver_explicit

  real(dp), parameter :: eps_H = eps_sp_dp
                         ! If computations are in double precision,
                         ! single precision is a good choice
                         ! for the ice-thickness epsilon

#if !defined(ALLOW_OPENAD) /* Normal */
  private
#endif
  public :: calc_thk_init
  public :: calc_thk_sia_expl, calc_thk_sia_impl
  public :: calc_thk_expl, calc_thk_impl
  public :: calc_thk_mask_update
  public :: account_mb_source
contains

!-------------------------------------------------------------------------------
!> Initialisations for the ice thickness computation.
!<------------------------------------------------------------------------------
subroutine calc_thk_init()

implicit none

#if defined(ALLOW_OPENAD) /* OpenAD */
integer(i4b) :: i, j
#endif /* OpenAD */

#if (defined(CALCZS))
  errormsg = ' >>> calc_thk_init: Replace CALCZS by CALCTHK in header file!'
  call error(errormsg)
#elif (!defined(CALCTHK))
  errormsg = ' >>> calc_thk_init: Define CALCTHK in header file!'
  call error(errormsg)
#endif

!-------- Computation/initialisation of the ice base topography
!                                       and its time derivative --------

#if (MARGIN==1 || MARGIN==2)   /* only grounded ice */

zb_new   = zl_new
dzb_dtau = dzl_dtau

#elif (MARGIN==3)   /* grounded and floating ice */

#if !defined(ALLOW_OPENAD) /* Normal */

where (mask <= 1_i1b)   ! grounded ice or ice-free land
   zb_new   = zl_new
   dzb_dtau = dzl_dtau
elsewhere   ! (mask >= 2_i1b; ocean or floating ice)
   zb_new   = zb       ! initialisation,
   dzb_dtau = 0.0_dp   ! will be overwritten later
end where

#else /* OpenAD */

do i=0, IMAX
do j=0, JMAX
   if (mask(j,i) <= 1_i1b) then    ! grounded ice or ice-free land
      zb_new(j,i)   = zl_new(j,i)
      dzb_dtau(j,i) = dzl_dtau(j,i)
   else
      zb_new(j,i)   = zb(j,i)       ! initialisation,
      dzb_dtau(j,i) = 0.0_dp        ! will be overwritten later
   endif
end do
end do

#endif /* Normal vs. OpenAD */

#else

errormsg = ' >>> calc_thk_init: MARGIN must be either 1, 2 or 3!'
call error(errormsg)

#endif

!-------- Initialisation of the ice thickness
!                               and surface topography --------

zs_new = zs   ! initialisation,
H_new  = H    ! will be overwritten later

!-------- Solver type --------

#if (CALCTHK==1 || CALCTHK==4)   /* explicit solver */

  flag_solver_explicit = .true.

#elif (CALCTHK==2 || CALCTHK==3 \
       || CALCTHK==5 || CALCTHK==6)   /* implicit solver */

  flag_solver_explicit = .false.

#else

  errormsg = ' >>> calc_thk_init: CALCTHK must be between 1 and 6!'
  call error(errormsg)

#endif

!-------- Source term for the ice thickness equation --------

mb_source = as_perp - Q_b_tot - calving

end subroutine calc_thk_init

!-------------------------------------------------------------------------------
!> Explicit solver for the diffusive SIA ice surface equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_sia_expl(time, dtime, dxi, deta, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar

integer(i4b)                       :: i, j
real(dp)                           :: azs2, azs3
real(dp), dimension(0:JMAX,0:IMAX) :: czs2, czs3

!-------- Abbreviations --------

azs2 = dtime/(dxi*dxi)
azs3 = dtime/(deta*deta)

czs2 = 0.0_dp
czs3 = 0.0_dp

do i=0, IMAX-1
do j=0, JMAX
   czs2(j,i) = azs2*0.5_dp*(h_diff(j,i)+h_diff(j,i+1)) &
               *(sq_g22_sgx(j,i)*insq_g11_sgx(j,i))
end do
end do

do i=0, IMAX
do j=0, JMAX-1
   czs3(j,i) = azs3*0.5_dp*(h_diff(j,i)+h_diff(j+1,i)) &
               *(sq_g11_sgy(j,i)*insq_g22_sgy(j,i))
end do
end do

!-------- Solution of the explicit scheme --------

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then   ! inner point

      zs_new(j,i) = zs(j,i) &
                    + dtime*(mb_source(j,i)+dzb_dtau(j,i)) &
                    + ( ( czs2(j,i)  *(zs(j,i+1)-zs(j,i)  ) &
                         -czs2(j,i-1)*(zs(j,i)  -zs(j,i-1)) ) &
                       +( czs3(j,i)  *(zs(j+1,i)-zs(j,i)  ) &
                         -czs3(j-1,i)*(zs(j,i)  -zs(j-1,i)) ) ) &
                      *(insq_g11_g(j,i)*insq_g22_g(j,i))

   else
      zs_new(j,i) = zb_new(j,i)   ! zero-thickness boundary condition
   end if

end do
end do

!-------- Ice thickness --------

H_new = zs_new - zb_new

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call thk_adjust_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_sia_expl

!-------------------------------------------------------------------------------
!> Over-implicit solver for the diffusive SIA ice surface equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_sia_impl(time, dtime, dxi, deta, z_mar, mean_accum)

#if !defined(ALLOW_OPENAD) /* Normal */
use sico_maths_m, only : sor_sprs
#else /* OpenAD */
use sico_maths_m
#endif /* Normal vs. OpenAD */

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar
real(dp), intent(in) :: mean_accum

integer(i4b)                       :: i, j
integer(i4b)                       :: k, nnz
real(dp)                           :: azs2, azs3
real(dp), dimension(0:JMAX,0:IMAX) :: czs2, czs3

#if (CALCTHK==2)

integer(i4b)                            :: ierr
integer(i4b)                            :: iter
integer(i4b)                            :: nc, nr
integer(i4b), parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)

#if !defined(ALLOW_OPENAD) /* Normal */
integer(i4b), allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
integer(i4b), allocatable, dimension(:) :: lgs_a_diag_index
integer(i4b), allocatable, dimension(:) :: lgs_a_index_trim
real(dp),     allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
real(dp),     allocatable, dimension(:) :: lgs_a_value_trim
#else /* OpenAD */
real(dp),             dimension(n_sprs) :: lgs_a_index
integer(i4b),         dimension(nmax+1) :: lgs_a_ptr
integer(i4b),           dimension(nmax) :: lgs_a_diag_index
integer(i4b),         dimension(n_sprs) :: lgs_a_index_trim
real(dp),             dimension(n_sprs) :: lgs_a_value
real(dp),               dimension(nmax) :: lgs_b_value
real(dp),               dimension(nmax) :: lgs_x_value
real(dp),             dimension(n_sprs) :: lgs_a_value_trim
#endif /* Normal vs. OpenAD */

real(dp)                                :: eps_sor

#elif (CALCTHK==3)

#if !defined(ALLOW_OPENAD) /* Normal */
LIS_INTEGER                            :: ierr
LIS_INTEGER                            :: iter
LIS_INTEGER                            :: nc, nr
LIS_INTEGER, parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
LIS_INTEGER, parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_diag_index
LIS_MATRIX                             :: lgs_a
LIS_VECTOR                             :: lgs_b, lgs_x
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
LIS_SOLVER                             :: solver
character(len=256)                     :: ch_solver_set_option
#else /* OpenAD */
integer(i4b)                            :: nc, nr
integer(i4b), parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)
real(dp),             dimension(n_sprs) :: lgs_a_index
integer(i4b),         dimension(nmax+1) :: lgs_a_ptr
integer(i4b),           dimension(nmax) :: lgs_a_diag_index
integer(i4b),         dimension(n_sprs) :: lgs_a_index_trim
real(dp),             dimension(n_sprs) :: lgs_a_value
real(dp),               dimension(nmax) :: lgs_b_value
real(dp),               dimension(nmax) :: lgs_x_value
real(dp),             dimension(n_sprs) :: lgs_a_value_trim
#endif /* Normal vs. OpenAD */

#endif

!-------- Abbreviations --------

azs2 = dtime/(dxi*dxi)
azs3 = dtime/(deta*deta)

czs2 = 0.0_dp
czs3 = 0.0_dp

do i=0, IMAX-1
do j=0, JMAX
  czs2(j,i) = azs2*0.5_dp*(h_diff(j,i)+h_diff(j,i+1)) &
              *(sq_g22_sgx(j,i)*insq_g11_sgx(j,i))
end do
end do

do i=0, IMAX
do j=0, JMAX-1
   czs3(j,i) = azs3*0.5_dp*(h_diff(j,i)+h_diff(j+1,i)) &
               *(sq_g11_sgy(j,i)*insq_g22_sgy(j,i))
end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

#if (CALCTHK==2 || CALCTHK==3)

#if !defined(ALLOW_OPENAD) /* Normal */
allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))
#endif /* Normal */

lgs_a_value = 0.0_dp
#if !defined(ALLOW_OPENAD) /* Normal */
lgs_a_index = 0
#else /* OpenAD */
lgs_a_index = 0.0_dp
#endif /* Normal vs. OpenAD */
lgs_a_ptr   = 0

lgs_b_value = 0.0_dp
lgs_x_value = 0.0_dp

lgs_a_ptr(1) = 1

k = 0

do nr=1, nmax   ! loop over rows

   i = n2i(nr)
   j = n2j(nr)

   if (flag_inner_point(j,i)) then   ! inner point

      nc = ij2n(j,i-1)   ! smallest nc (column counter), for zs(j,i-1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs2(j,i-1)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      nc = ij2n(j-1,i)   ! next nc (column counter), for zs(j-1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs3(j-1,i)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      nc = ij2n(j,i)     ! next nc (column counter), for zs(j,i)
      ! if (nc /= nr) &                     (diagonal element)
      !    stop ' >>> calc_thk_sia_impl: Check for diagonal element failed!'
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = 1.0_dp &
                       + ((czs2(j,i)+czs2(j,i-1))+(czs3(j,i)+czs3(j-1,i))) &
                         *OVI_WEIGHT &
                         *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_diag_index(nr) = k
      lgs_a_index(k) = nc

      nc = ij2n(j+1,i)   ! next nc (column counter), for zs(j+1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs3(j,i)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      nc = ij2n(j,i+1)   ! largest nc (column counter), for zs(j,i+1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs2(j,i)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      lgs_b_value(nr) = zs(j,i) &
                          + dtime*(mb_source(j,i)+dzb_dtau(j,i)) &
                          + ( ( czs2(j,i)*(zs(j,i+1)-zs(j,i)) &
                               -czs2(j,i-1)*(zs(j,i)-zs(j,i-1)) ) &
                             +( czs3(j,i)*(zs(j+1,i)-zs(j,i)) &
                               -czs3(j-1,i)*(zs(j,i)-zs(j-1,i)) ) ) &
                            *(1.0_dp-OVI_WEIGHT) &
                            *(insq_g11_g(j,i)*insq_g22_g(j,i))
                                                          ! right-hand side

   else   !  zero-thickness boundary condition

      k = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k)       = 1.0_dp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = zb_new(j,i)

   end if

   lgs_x_value(nr) = zs(j,i)   ! old surface topography,
                               ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

#else

errormsg = ' >>> calc_thk_sia_impl: CALCTHK must be either 2 or 3!'
call error(errormsg)

#endif

!-------- Solution of the system of linear equations --------

#if (CALCTHK==2)

!  ------ Solution with the built-in SOR solver

#if !defined(ALLOW_OPENAD) /* Normal */
allocate(lgs_a_value_trim(nnz), lgs_a_index_trim(nnz))
#endif /* Normal */

do k=1, nnz   ! relocate matrix to trimmed arrays
   lgs_a_value_trim(k) = lgs_a_value(k)
#if !defined(ALLOW_OPENAD) /* Normal */
   lgs_a_index_trim(k) = lgs_a_index(k)
#else /* OpenAD */
   lgs_a_index_trim(k) = int(lgs_a_index(k))
#endif /* Normal vs. OpenAD */
end do

#if !defined(ALLOW_OPENAD) /* Normal */
deallocate(lgs_a_value, lgs_a_index)
#endif /* Normal */

eps_sor = 1.0e-05_dp*mean_accum*dtime   ! convergence parameter

call sor_sprs(lgs_a_value_trim, &
              lgs_a_index_trim, lgs_a_diag_index, lgs_a_ptr, &
              lgs_b_value, &
              nnz, nmax, &
#if defined(ALLOW_OPENAD) /* OpenAD */
              n_sprs, &
#endif /* OpenAD */
              OMEGA_SOR, eps_sor, lgs_x_value, ierr)

do nr=1, nmax
   i = n2i(nr)
   j = n2j(nr)
   zs_new(j,i) = lgs_x_value(nr)
end do

#if !defined(ALLOW_OPENAD) /* Normal */
deallocate(lgs_a_value_trim, lgs_a_index_trim, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)
#endif /* Normal */

#elif (CALCTHK==3)

#if !defined(ALLOW_OPENAD) /* Normal */
  
  !  ------ Settings for Lis
  
  call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
  call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
  call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)
  
  call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
  call lis_vector_set_size(lgs_b, 0, nmax, ierr)
  call lis_vector_set_size(lgs_x, 0, nmax, ierr)
  
  do nr=1, nmax
  
     do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
        call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                                 lgs_a_value(nc), lgs_a, ierr)
     end do
  
     call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
     call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)
  
  end do

  call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
  call lis_matrix_assemble(lgs_a, ierr)
  
  !  ------ Solution with Lis
  
  call lis_solver_create(solver, ierr)
  
  ch_solver_set_option = '-i bicg -p ilu '// &
                         '-maxiter 1000 -tol 1.0e-12 -initx_zeros false'
  
  call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
  call CHKERR(ierr)
  
  call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
  call CHKERR(ierr)
  
  call lis_solver_get_iter(solver, iter, ierr)
  write(6,'(10x,a,i0)') 'calc_thk_sia_impl: iter = ', iter
  
  lgs_x_value = 0.0_dp
  call lis_vector_gather(lgs_x, lgs_x_value, ierr)
  call lis_matrix_destroy(lgs_a, ierr)
  call lis_vector_destroy(lgs_b, ierr)
  call lis_vector_destroy(lgs_x, ierr)
  call lis_solver_destroy(solver, ierr)

  do nr=1, nmax
    i = n2i(nr)
    j = n2j(nr)
    zs_new(j,i) = lgs_x_value(nr)
  end do
  
  deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
  deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)
  
#else /* OpenAD */

  do k=1, nnz   ! relocate matrix to trimmed arrays
     lgs_a_value_trim(k) = lgs_a_value(k)
     lgs_a_index_trim(k) = int(lgs_a_index(k))
  end do

  call sico_lis_solver(nmax, nnz, &
                       lgs_a_ptr, lgs_a_index_trim, &
                       lgs_a_value_trim, lgs_b_value, lgs_x_value)

  do nr=1, nmax
     i = n2i(nr)
     j = n2j(nr)
     zs_new(j,i) = lgs_x_value(nr)
  end do

#endif /* Normal vs. OpenAD */

#endif

!-------- Ice thickness --------

H_new = zs_new - zb_new

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call thk_adjust_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_sia_impl

!-------------------------------------------------------------------------------
!> Explicit solver for the general ice thickness equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_expl(time, dtime, dxi, deta, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar

integer(i4b)                       :: i, j
real(dp), dimension(0:JMAX,0:IMAX) :: dt_darea
real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_1, vx_m_2, vy_m_1, vy_m_2
real(dp), dimension(0:JMAX,0:IMAX) :: upH_x_1, upH_x_2, upH_y_1, upH_y_2
real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_x_1, sq_g22_x_2
real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_y_1, sq_g11_y_2

!-------- Abbreviations --------

dt_darea = dtime/area

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then

      vx_m_1(j,i) = vx_m(j,i-1)
      vx_m_2(j,i) = vx_m(j,i)
      vy_m_1(j,i) = vy_m(j-1,i)
      vy_m_2(j,i) = vy_m(j,i)

      if (vx_m_1(j,i) >= 0.0_dp) then
         upH_x_1(j,i) = H(j,i-1)
      else
         upH_x_1(j,i) = H(j,i)
      end if

      if (vx_m_2(j,i) >= 0.0_dp) then
         upH_x_2(j,i) = H(j,i)
      else
         upH_x_2(j,i) = H(j,i+1)
      end if

      if (vy_m_1(j,i) >= 0.0_dp) then
         upH_y_1(j,i) = H(j-1,i)
      else
         upH_y_1(j,i) = H(j,i)
      end if

      if (vy_m_2(j,i) >= 0.0_dp) then
         upH_y_2(j,i) = H(j,i)
      else
         upH_y_2(j,i) = H(j+1,i)
      end if

      sq_g22_x_1(j,i) = sq_g22_sgx(j,i-1)
      sq_g22_x_2(j,i) = sq_g22_sgx(j,i)
      sq_g11_y_1(j,i) = sq_g11_sgy(j-1,i)
      sq_g11_y_2(j,i) = sq_g11_sgy(j,i)

   else   ! .not.(flag_inner_point(j,i))

      vx_m_1(j,i) = 0.0_dp
      vx_m_2(j,i) = 0.0_dp
      vy_m_1(j,i) = 0.0_dp
      vy_m_2(j,i) = 0.0_dp

      upH_x_1(j,i) = 0.0_dp
      upH_x_2(j,i) = 0.0_dp
      upH_y_1(j,i) = 0.0_dp
      upH_y_2(j,i) = 0.0_dp

      sq_g22_x_1(j,i) = 0.0_dp
      sq_g22_x_2(j,i) = 0.0_dp
      sq_g11_y_1(j,i) = 0.0_dp
      sq_g11_y_2(j,i) = 0.0_dp

   end if

end do
end do

!-------- Solution of the explicit scheme --------

#if !defined(ALLOW_OPENAD) /* Normal */

where (flag_inner_point)   ! inner point

   H_new = H + dtime*mb_source &
             - dt_darea &
               * (  ( vx_m_2*upH_x_2*sq_g22_x_2*deta   &
                     -vx_m_1*upH_x_1*sq_g22_x_1*deta ) &
                  + ( vy_m_2*upH_y_2*sq_g11_y_2*dxi    &
                     -vy_m_1*upH_y_1*sq_g11_y_1*dxi  ) )

elsewhere
   H_new = 0.0_dp   ! zero-thickness boundary condition
end where

#else /* OpenAD */

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then   ! inner point

      H_new(j,i) = H(j,i) + dtime*mb_source(j,i) &
                - dt_darea(j,i) &
                  * (  ( vx_m_2(j,i)*upH_x_2(j,i)*sq_g22_x_2(j,i)*deta   &
                        -vx_m_1(j,i)*upH_x_1(j,i)*sq_g22_x_1(j,i)*deta ) &
                     + ( vy_m_2(j,i)*upH_y_2(j,i)*sq_g11_y_2(j,i)*dxi    &
                        -vy_m_1(j,i)*upH_y_1(j,i)*sq_g11_y_1(j,i)*dxi  ) )

   else
      H_new(j,i) = 0.0_dp   ! zero-thickness boundary condition
   end if

end do
end do

#endif /* Normal vs. OpenAD */

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call thk_adjust_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_expl

!-------------------------------------------------------------------------------
!> Over-implicit solver for the general ice thickness equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_impl(time, dtime, dxi, deta, z_mar, mean_accum)

#if !defined(ALLOW_OPENAD) /* Normal */
#if (CALCTHK==5)
use sico_maths_m, only : sor_sprs
#elif (CALCTHK==6)
use sico_maths_m, only : sor_sprs, sico_lis_solver
#endif
#else /* OpenAD */
use sico_maths_m
#endif /* Normal vs. OpenAD */

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar
real(dp), intent(in) :: mean_accum

integer(i4b)                       :: i, j
integer(i4b)                       :: k, nnz
real(dp), dimension(0:JMAX,0:IMAX) :: dt_darea
real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_1, vx_m_2, vy_m_1, vy_m_2
real(dp), dimension(0:JMAX,0:IMAX) :: upH_x_1, upH_x_2, upH_y_1, upH_y_2
real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_x_1, sq_g22_x_2
real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_y_1, sq_g11_y_2

#if (CALCTHK==5)

integer(i4b)                            :: ierr
integer(i4b)                            :: iter
integer(i4b)                            :: nc, nr
integer(i4b), parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)

#if !defined(ALLOW_OPENAD) /* Normal */
integer(i4b), allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
integer(i4b), allocatable, dimension(:) :: lgs_a_diag_index
integer(i4b), allocatable, dimension(:) :: lgs_a_index_trim
real(dp),     allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
real(dp),     allocatable, dimension(:) :: lgs_a_value_trim
#else /* OpenAD */
real(dp),             dimension(n_sprs) :: lgs_a_index
integer(i4b),         dimension(nmax+1) :: lgs_a_ptr
integer(i4b),           dimension(nmax) :: lgs_a_diag_index
integer(i4b),         dimension(n_sprs) :: lgs_a_index_trim
real(dp),             dimension(n_sprs) :: lgs_a_value
real(dp),               dimension(nmax) :: lgs_b_value
real(dp),               dimension(nmax) :: lgs_x_value
real(dp),             dimension(n_sprs) :: lgs_a_value_trim
#endif /* Normal vs. OpenAD */

real(dp)                                :: eps_sor

#elif (CALCTHK==6)

#if !defined(ALLOW_OPENAD) /* Normal */
LIS_INTEGER                            :: ierr
LIS_INTEGER                            :: iter
LIS_INTEGER                            :: nc, nr
LIS_INTEGER, parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
LIS_INTEGER, parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_diag_index
LIS_MATRIX                             :: lgs_a
LIS_VECTOR                             :: lgs_b, lgs_x
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
LIS_SOLVER                             :: solver
character(len=256)                     :: ch_solver_set_option
#else /* OpenAD */
integer(i4b)                            :: nc, nr
integer(i4b), parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)
real(dp),             dimension(n_sprs) :: lgs_a_index
integer(i4b),         dimension(nmax+1) :: lgs_a_ptr
integer(i4b),           dimension(nmax) :: lgs_a_diag_index
integer(i4b),         dimension(n_sprs) :: lgs_a_index_trim
real(dp),             dimension(n_sprs) :: lgs_a_value
real(dp),               dimension(nmax) :: lgs_b_value
real(dp),               dimension(nmax) :: lgs_x_value
real(dp),             dimension(n_sprs) :: lgs_a_value_trim
#endif /* Normal vs. OpenAD */

#endif

!-------- Abbreviations --------

dt_darea = dtime/area

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then

      vx_m_1(j,i) = vx_m(j,i-1)
      vx_m_2(j,i) = vx_m(j,i)
      vy_m_1(j,i) = vy_m(j-1,i)
      vy_m_2(j,i) = vy_m(j,i)

      if (vx_m_1(j,i) >= 0.0_dp) then
         upH_x_1(j,i) = H(j,i-1)
      else
         upH_x_1(j,i) = H(j,i)
      end if

      if (vx_m_2(j,i) >= 0.0_dp) then
         upH_x_2(j,i) = H(j,i)
      else
         upH_x_2(j,i) = H(j,i+1)
      end if

      if (vy_m_1(j,i) >= 0.0_dp) then
         upH_y_1(j,i) = H(j-1,i)
      else
         upH_y_1(j,i) = H(j,i)
      end if

      if (vy_m_2(j,i) >= 0.0_dp) then
         upH_y_2(j,i) = H(j,i)
      else
         upH_y_2(j,i) = H(j+1,i)
      end if

      sq_g22_x_1(j,i) = sq_g22_sgx(j,i-1)
      sq_g22_x_2(j,i) = sq_g22_sgx(j,i)
      sq_g11_y_1(j,i) = sq_g11_sgy(j-1,i)
      sq_g11_y_2(j,i) = sq_g11_sgy(j,i)

   else   ! .not.(flag_inner_point(j,i))

      vx_m_1(j,i) = 0.0_dp
      vx_m_2(j,i) = 0.0_dp
      vy_m_1(j,i) = 0.0_dp
      vy_m_2(j,i) = 0.0_dp

      upH_x_1(j,i) = 0.0_dp
      upH_x_2(j,i) = 0.0_dp
      upH_y_1(j,i) = 0.0_dp
      upH_y_2(j,i) = 0.0_dp

      sq_g22_x_1(j,i) = 0.0_dp
      sq_g22_x_2(j,i) = 0.0_dp
      sq_g11_y_1(j,i) = 0.0_dp
      sq_g11_y_2(j,i) = 0.0_dp

   end if

end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

#if (CALCTHK==5 || CALCTHK==6)

#if !defined(ALLOW_OPENAD) /* Normal */
allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))
#endif /* Normal */

lgs_a_value = 0.0_dp
#if !defined(ALLOW_OPENAD) /* Normal */
lgs_a_index = 0
#else /* OpenAD */
lgs_a_index = 0.0_dp
#endif /* Normal vs. OpenAD */
lgs_a_ptr   = 0

lgs_b_value = 0.0_dp
lgs_x_value = 0.0_dp

lgs_a_ptr(1) = 1

k = 0

do nr=1, nmax   ! loop over rows

   i = n2i(nr)
   j = n2j(nr)

   if (flag_inner_point(j,i)) then   ! inner point

      k=k+1 ; nc=ij2n(j,i-1) ; lgs_a_index(k)=nc   ! for H(j,i-1)
      if (vx_m_1(j,i) > 0.0_dp) &
         lgs_a_value(k) = -dt_darea(j,i)*vx_m_1(j,i) &
                                        *sq_g22_x_1(j,i)*deta*OVI_WEIGHT

      k=k+1 ; nc=ij2n(j-1,i) ; lgs_a_index(k)=nc   ! for H(j-1,i)
      if (vy_m_1(j,i) > 0.0_dp) &
         lgs_a_value(k) = -dt_darea(j,i)*vy_m_1(j,i) &
                                        *sq_g11_y_1(j,i)*dxi*OVI_WEIGHT

      k=k+1 ; lgs_a_index(k)=nr ; lgs_a_diag_index(nr)=k  ! for H(j,i)
      lgs_a_value(k) = 1.0_dp                             ! (diagonal element)
      if (vy_m_1(j,i) < 0.0_dp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          - dt_darea(j,i)*vy_m_1(j,i) &
                                         *sq_g11_y_1(j,i)*dxi*OVI_WEIGHT
      if (vx_m_1(j,i) < 0.0_dp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          - dt_darea(j,i)*vx_m_1(j,i) &
                                         *sq_g22_x_1(j,i)*deta*OVI_WEIGHT
      if (vx_m_2(j,i) > 0.0_dp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          + dt_darea(j,i)*vx_m_2(j,i) &
                                         *sq_g22_x_2(j,i)*deta*OVI_WEIGHT
      if (vy_m_2(j,i) > 0.0_dp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          + dt_darea(j,i)*vy_m_2(j,i) &
                                         *sq_g11_y_2(j,i)*dxi*OVI_WEIGHT

      k=k+1 ; nc=ij2n(j+1,i) ; lgs_a_index(k)=nc   ! for H(j+1,i)
      if (vy_m_2(j,i) < 0.0_dp) &
         lgs_a_value(k) = dt_darea(j,i)*vy_m_2(j,i) &
                                       *sq_g11_y_2(j,i)*dxi*OVI_WEIGHT

      k=k+1 ; nc=ij2n(j,i+1) ; lgs_a_index(k)=nc   ! for H(j,i+1)
      if (vx_m_2(j,i) < 0.0_dp) &
         lgs_a_value(k) = dt_darea(j,i)*vx_m_2(j,i) &
                                       *sq_g22_x_2(j,i)*deta*OVI_WEIGHT

      lgs_b_value(nr) = H(j,i) &
                        +dtime*mb_source(j,i) &
                        -(1.0_dp-OVI_WEIGHT) &
                           * dt_darea(j,i) &
                             * (  ( vx_m_2(j,i)*upH_x_2(j,i) &
                                               *sq_g22_x_2(j,i)*deta   &
                                   -vx_m_1(j,i)*upH_x_1(j,i) &
                                               *sq_g22_x_1(j,i)*deta ) &
                                + ( vy_m_2(j,i)*upH_y_2(j,i) &
                                               *sq_g11_y_2(j,i)*dxi    &
                                   -vy_m_1(j,i)*upH_y_1(j,i) &
                                               *sq_g11_y_1(j,i)*dxi  ) )
                                                          ! right-hand side

   else   ! zero-thickness boundary condition

      k = k+1
      lgs_a_value(k)       = 1.0_dp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = 0.0_dp

   end if

   lgs_x_value(nr) = H(j,i)   ! old ice thickness,
                              ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

#else

errormsg = ' >>> calc_thk_impl: CALCTHK must be either 5 or 6!'
call error(errormsg)

#endif

!-------- Solution of the system of linear equations --------

#if (CALCTHK==5)

!  ------ Solution with the built-in SOR solver

#if !defined(ALLOW_OPENAD) /* Normal */
allocate(lgs_a_value_trim(nnz), lgs_a_index_trim(nnz))
#endif /* Normal */

do k=1, nnz   ! relocate matrix to trimmed arrays
   lgs_a_value_trim(k) = lgs_a_value(k)
#if !defined(ALLOW_OPENAD) /* Normal */
   lgs_a_index_trim(k) = lgs_a_index(k)
#else /* OpenAD */
   lgs_a_index_trim(k) = int(lgs_a_index(k))
#endif /* Normal vs. OpenAD */
end do

#if !defined(ALLOW_OPENAD) /* Normal */
deallocate(lgs_a_value, lgs_a_index)
#endif /* Normal */

eps_sor = 1.0e-05_dp*mean_accum*dtime   ! convergence parameter

call sor_sprs(lgs_a_value_trim, &
              lgs_a_index_trim, lgs_a_diag_index, lgs_a_ptr, &
              lgs_b_value, &
              nnz, nmax, &
#if defined(ALLOW_OPENAD) /* OpenAD */
              n_sprs, &
#endif /* OpenAD */
              OMEGA_SOR, eps_sor, lgs_x_value, ierr)

do nr=1, nmax
   i = n2i(nr)
   j = n2j(nr)
   H_new(j,i) = lgs_x_value(nr)
end do

#if !defined(ALLOW_OPENAD) /* Normal */
deallocate(lgs_a_value_trim, lgs_a_index_trim, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)
#endif /* Normal */

#elif (CALCTHK==6)

#if !defined(ALLOW_OPENAD) /* Normal */

!  ------ Settings for Lis

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
call lis_matrix_assemble(lgs_a, ierr)

!  ------ Solution with Lis

call lis_solver_create(solver, ierr)

ch_solver_set_option = '-i bicg -p ilu '// &
                       '-maxiter 1000 -tol 1.0e-12 -initx_zeros false'

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

call lis_solver_get_iter(solver, iter, ierr)
write(6,'(10x,a,i0)') 'calc_thk_impl: iter = ', iter

lgs_x_value = 0.0_dp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

do nr=1, nmax
   i = n2i(nr)
   j = n2j(nr)
   H_new(j,i) = lgs_x_value(nr)
end do

deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)

#else /* OpenAD */

do k=1, nnz   ! relocate matrix to trimmed arrays
   lgs_a_value_trim(k) = lgs_a_value(k)
#if !defined(ALLOW_OPENAD) /* Normal */
   lgs_a_index_trim(k) = lgs_a_index(k)
#else /* OpenAD */
   lgs_a_index_trim(k) = int(lgs_a_index(k))
#endif /* Normal vs. OpenAD */
end do
call sico_lis_solver(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index_trim, &
                           lgs_a_value_trim, lgs_b_value, lgs_x_value)
do nr=1, nmax
   i = n2i(nr)
   j = n2j(nr)
   H_new(j,i) = lgs_x_value(nr)
end do

#endif /* Normal vs. OpenAD */

#endif

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call thk_adjust_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_impl

!-------------------------------------------------------------------------------
!> Ice thickness evolution due to the source term (surface mass balance,
!! basal mass balance and calving).
!<------------------------------------------------------------------------------
subroutine apply_mb_source(dtime, z_mar)

  implicit none

  real(dp), intent(in) :: dtime
  real(dp), intent(in) :: z_mar

  integer(i4b) :: i, j

!-------- Compute new ice thickness H_new_flow due to glacial flow only
!         (no source term considered) --------

  H_new_flow = 0.0_dp

  do i=0, IMAX
  do j=0, JMAX
     if (flag_inner_point(j,i)) then
        H_new_flow(j,i) = H_new(j,i) - dtime*mb_source(j,i)
                          ! new ice thickness due to glacial flow only
     end if
  end do
  end do

!-------- Apply source term --------

  H_new = 0.0_dp

#if (MB_ACCOUNT==0)

  do i=0, IMAX
  do j=0, JMAX
     if (flag_inner_point(j,i)) then
        H_new(j,i) = max(H_new_flow(j,i) + dtime*mb_source(j,i), 0.0_dp)
     end if
  end do
  end do

#elif (MB_ACCOUNT==1)

  do i=2, IMAX-2   ! outermost two grid points excepted,
  do j=2, JMAX-2   ! required for accurate mass balance accounting
     H_new(j,i) = max(H_new_flow(j,i) + dtime*mb_source(j,i), 0.0_dp)
  end do
  end do

#else

  errormsg = ' >>> apply_mb_source: MB_ACCOUNT must be either 0 or 1!'
  call error(errormsg)

#endif

!-------- Reset new ice thickness if needed --------

#if (MARGIN==1)

  do i=0, IMAX
  do j=0, JMAX
     if (mask(j,i) >= 2_i1b .and. H_new(j,i) > eps_H) then
        H_new(j,i) = 0.0_dp
     end if
  end do
  end do

#else   /* MARGIN==2 or 3  */

  do i=0, IMAX
  do j=0, JMAX
     if (zl(j,i) <= Z_ABYSS) then
        H_new(j,i) = 0.0_dp
     end if
  end do
  end do

#endif

end subroutine apply_mb_source

!-------------------------------------------------------------------------------
!> Adjustment of the newly computed ice thickness distribution.
!<------------------------------------------------------------------------------
  subroutine thk_adjust(time, dtime)

  implicit none

  real(dp), intent(in) :: time, dtime

  real(dp), dimension(0:JMAX,0:IMAX) :: H_new_tmp
  real(dp)                           :: time_dimless
  real(dp)                           :: target_topo_tau_factor, target_topo_tau
  real(dp)                           :: dtime_inv

#if defined(ALLOW_OPENAD) /* OpenAD */
  integer(i4b) :: i, j
#endif /* OpenAD */

!-------- Saving computed H_new before any adjustments --------

  H_new_tmp = H_new

!-------- Correct negative thickness values --------

#if !defined(ALLOW_OPENAD) /* Normal */

  where (H_new < 0.0_dp) H_new = 0.0_dp

#else /* OpenAD */

  do i=0, IMAX  
  do j=0, JMAX 
     if ( H_new(j,i) < 0.0_dp ) then
        H_new(j,i) = 0.0_dp
     end if
  end do
  end do

#endif /* Normal vs. OpenAD */

!-------- Further adjustments --------

#if (THK_EVOL==0)

!  ------ No evolution of the ice thickness

  H_new = H   ! newly computed ice thickness is discarded

#elif (THK_EVOL==1)

!  ------ Ice thickness evolves freely

  !!! continue

#elif (THK_EVOL==2)

!  ------ Nudging towards prescribed target topography
!                                    with varying relaxation time

  if (time >= time_target_topo_final) then

     H_new = H_target

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

     H_new =   ( target_topo_tau*H_new + dtime*H_target ) &
             / ( target_topo_tau       + dtime )

  end if

#elif (THK_EVOL==3)

!  ------ Nudging towards prescribed target topography
!                                    with constant relaxation time

  target_topo_tau = target_topo_tau_0

  H_new =   ( target_topo_tau*H_new + dtime*H_target ) &
          / ( target_topo_tau       + dtime )

#elif (THK_EVOL==4)

!  ------ Maximum ice extent constrained by prescribed mask

#if !defined(ALLOW_OPENAD) /* Normal */

  where (mask_maxextent == 0_i1b) &   ! not allowed to glaciate
     H_new = 0.0_dp

#else /* OpenAD */

  do i=0, IMAX
  do j=0, JMAX
    if ( mask_maxextent(j,i) == 0_i1b ) then
      H_new(j,i) = 0.0_dp
    end if
  end do
  end do

#endif /* Normal vs. OpenAD */

#else

  errormsg = ' >>> thk_adjust: THK_EVOL must be between 0 and 4!'
  call error(errormsg)

#endif

!-------- Computation of the mass balance adjustment --------

  dtime_inv = 1.0_dp/dtime

  smb_corr = (H_new-H_new_tmp)*dtime_inv

  as_perp = as_perp + smb_corr
  accum   = accum   + max(smb_corr, 0.0_dp)
  runoff  = runoff  - min(smb_corr, 0.0_dp)
                    ! runoff is counted as positive for mass loss

end subroutine thk_adjust

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
!-------------------------------------------------------------------------------
!> Adjustment of the newly computed ice thickness distribution due to either
!  the retreat mask due to oceanic forcing or the ice-shelf collapse mask
!  (counted as calving).
!<------------------------------------------------------------------------------
  subroutine thk_adjust_retreat_mask(time, dtime)

  implicit none

  real(dp), intent(in) :: time, dtime

  integer(i4b) :: i, j
  real(dp), dimension(0:JMAX,0:IMAX) :: H_new_tmp, dHdt_retreat
  real(dp)                           :: dtime_inv
  real(dp)                           :: dtime_1year, dtime_1year_inv

  dtime_inv       = 1.0_dp/dtime
  dtime_1year     = year2sec   ! 1 year (in seconds)
  dtime_1year_inv = 1.0_dp/dtime_1year

!-------- Saving computed H_new before any adjustments --------

  H_new_tmp = H_new

!-------- Adjustment due to the retreat mask --------

  dHdt_retreat = 0.0_dp

  do i=0, IMAX
  do j=0, JMAX

#if (RETREAT_MASK==1)
     if (H_new(j,i) > 0.0_dp) then
#elif (ICE_SHELF_COLLAPSE_MASK==1)
     if ((H_new(j,i) > 0.0_dp).and.(mask(j,i)==3_i1b)) then
#endif

        dHdt_retreat(j,i) = -(1.0_dp-r_mask_retreat(j,i))*H_ref_retreat(j,i) &
                                                         *dtime_1year_inv

        H_new(j,i) = max((H_new(j,i) + dHdt_retreat(j,i)*dtime), 0.0_dp)

     end if

  end do
  end do

!-------- Computation of the mass balance adjustment --------

  smb_corr_retreat_mask = (H_new-H_new_tmp)*dtime_inv

  calving = calving - smb_corr_retreat_mask
                    ! calving is counted as positive for mass loss

end subroutine thk_adjust_retreat_mask

#endif   /* (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1) */

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask etc.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update(time, dtime, dxi, deta, z_sl, z_mar, &
                                n_calc_thk_mask_update_aux)

use topograd_m

implicit none

real(dp),     intent(in) :: time
real(dp),     intent(in) :: dtime, dxi, deta
real(dp),     intent(in) :: z_sl, z_mar
integer(i1b), intent(in) :: n_calc_thk_mask_update_aux

integer(i4b)                       :: i, j
real(dp)                           :: year_sec_inv
real(dp), dimension(0:JMAX,0:IMAX) :: H_new_tmp
real(dp)                           :: dtime_inv

!-------- Term abbreviations --------

year_sec_inv = 1.0_dp/year2sec

!-------- Saving computed H_new before any modifications --------

H_new_tmp = H_new

!-------- Update of the mask --------

dtime_inv = 1.0_dp/dtime

if (n_calc_thk_mask_update_aux == 1_i1b) then
   call calc_thk_mask_update_aux1(time, dtime, dtime_inv, z_sl, z_mar)
else if (n_calc_thk_mask_update_aux == 2_i1b) then
   call calc_thk_mask_update_aux2(time, dtime, dtime_inv, z_sl, z_mar)
else if (n_calc_thk_mask_update_aux == 3_i1b) then
   call calc_thk_mask_update_aux3(time, dtime, dtime_inv, z_sl, z_mar)
else
   errormsg = ' >>> calc_thk_mask_update:' &
            //         end_of_line &
            //'        n_calc_thk_mask_update_aux has no valid value!'
   call error(errormsg)
end if

!  ------ Enforce connectivity of the ocean

call ocean_connect()

!-------- Time derivatives --------

dzs_dtau  = (zs_new-zs)*dtime_inv
dzb_dtau  = (zb_new-zb)*dtime_inv
dzm_dtau  = dH_t_dtau+dzb_dtau
dH_dtau   = (H_new-H)*dtime_inv
dH_c_dtau = dzs_dtau-dzm_dtau

#if (THK_EVOL==2)
if ( abs((time-time_target_topo_final)*year_sec_inv) < eps ) then
   dzs_dtau  = 0.0_dp   ! Introduced
   dzb_dtau  = 0.0_dp   ! by
   dzm_dtau  = 0.0_dp   ! Tatsuru Sato
   dH_dtau   = 0.0_dp   ! for
   dH_c_dtau = 0.0_dp   ! stability
   dH_t_dtau = 0.0_dp   ! reasons
end if
#endif

!-------- New gradients --------

#if (TOPOGRAD==0)
call topograd_1(dxi, deta, 2)
#elif (TOPOGRAD==1)
call topograd_2(dxi, deta, 2)
#endif

!-------- Compute the volume flux across the CTS, am_perp --------

#if (CALCMOD==1)

do i=1, IMAX-1
do j=1, JMAX-1

   if ( ((mask(j,i) == 0_i1b).or.(mask(j,i) == 3_i1b)) &
        .and.(n_cts(j,i) == 1_i1b)) then

      am_perp_st(j,i) = &
              (   0.5_dp*(vx_c(0,j,i)+vx_c(0,j,i-1))*dzm_dxi_g(j,i) &
                + 0.5_dp*(vy_c(0,j,i)+vy_c(0,j-1,i))*dzm_deta_g(j,i) ) &
                - 0.5_dp*(vz_c(0,j,i)+vz_t(KTMAX-1,j,i))
      am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

   else
      am_perp_st(j,i) = 0.0_dp
      am_perp(j,i)    = 0.0_dp
   end if

end do
end do

#endif

end subroutine calc_thk_mask_update

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for SIA-only dynamics of ice sheets
!! without ice shelves.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux1(time, dtime, dtime_inv, z_sl, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dtime_inv, z_sl, z_mar

integer(i4b) :: i, j, kc, kt
real(dp)     :: rhosw_rho_ratio, rho_rhosw_ratio
real(dp)     :: H_inv

real(dp), dimension(0:JMAX,0:IMAX) :: H_sea_new, H_balance

!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW

!-------- Update of the mask --------

zs_new    = zb_new + H_new   ! ice surface topography
H_sea_new = z_sl - zl_new    ! sea depth
H_balance = H_sea_new*rhosw_rho_ratio

#if (THK_EVOL>=1)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1_i1b) then   ! grounded ice or ice-free land

      if (H_new(j,i) >= eps_H) then
         mask(j,i) = 0_i1b   ! grounded ice
      else
         mask(j,i) = 1_i1b   ! ice-free land
      end if

#if (MARGIN==2)

   else   ! (mask(j,i) == 2_i1b); sea

      if (H_new(j,i) >= eps_H) then

         if ( H_new(j,i) < (rhosw_rho_ratio*H_sea_new(j,i)) ) then

#if (MARINE_ICE_FORMATION==1)
            mask(j,i) = 2_i1b   ! floating ice cut off -> sea
#elif (MARINE_ICE_FORMATION==2)
            mask(j,i) = 0_i1b   ! "underwater ice"
#endif
         else
            mask(j,i) = 0_i1b   ! grounded ice
         end if

#if (MARINE_ICE_CALVING==2 \
      || MARINE_ICE_CALVING==4 \
      || MARINE_ICE_CALVING==6)
         if (zl0(j,i) < z_mar) mask(j,i) = 2_i1b   ! sea
#elif (MARINE_ICE_CALVING==3 \
        || MARINE_ICE_CALVING==5 \
        || MARINE_ICE_CALVING==7)
         if (zl_new(j,i) < z_mar) mask(j,i) = 2_i1b   ! sea
#endif

      else
         mask(j,i) = 2_i1b   ! sea
      end if

#endif

   end if

end do
end do

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2)
if (time >= time_target_topo_final) mask = mask_target
#endif

!-------- Correction of zs_new and zb_new for ice-free land and sea --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 1_i1b) then   ! ice-free land

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 2_i1b) then   ! sea

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 3_i1b) then   ! floating ice

      errormsg = ' >>> calc_thk_mask_update_aux1:' &
               //           end_of_line &
               //'          mask(j,i)==3 not allowed for' &
               //           end_of_line &
               //'          SIA-only dynamics!'
      call error(errormsg)

   end if

end do
end do

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice(z_sl)

!-------- Computation of further quantities --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 0_i1b) then   ! grounded ice

      if (n_cts(j,i) == 1_i1b) then
         if (H(j,i) > 0.0_dp) then
            H_inv        = 1.0_dp/H(j,i)
            H_c_new(j,i) = H_c(j,i) * H_new(j,i)*H_inv
            H_t_new(j,i) = H_t(j,i) * H_new(j,i)*H_inv
         else
            H_c_new(j,i) = 0.99_dp * H_new(j,i)   ! this case should not occur,
            H_t_new(j,i) = 0.01_dp * H_new(j,i)   ! just for safety
         end if
         zm_new(j,i) = zb_new(j,i)+H_t_new(j,i)
      else
         H_c_new(j,i) = H_new(j,i)
         H_t_new(j,i) = 0.0_dp
         zm_new(j,i)  = zb_new(j,i)
      end if

   else   ! mask(j,i) == 1_i1b or 2_i1b, ice-free land or sea

      H_c_new(j,i) = 0.0_dp
      H_t_new(j,i) = 0.0_dp
      H_new(j,i)   = 0.0_dp
      zs_new(j,i)  = zb_new(j,i)
      zm_new(j,i)  = zb_new(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux1

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for hybrid SIA/SStA dynamics of ice sheets
!! without ice shelves.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux2(time, dtime, dtime_inv, z_sl, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dtime_inv, z_sl, z_mar

integer(i4b) :: i, j, kc, kt
real(dp)     :: rhosw_rho_ratio, rho_rhosw_ratio
real(dp)     :: H_inv

real(dp), dimension(0:JMAX,0:IMAX) :: H_sea_new, H_balance

!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW

!-------- Update of the mask --------

zs_new    = zb_new + H_new   ! ice surface topography
H_sea_new = z_sl - zl_new    ! sea depth
H_balance = H_sea_new*rhosw_rho_ratio

#if (THK_EVOL>=1)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1_i1b) then   ! grounded ice or ice-free land

      if (H_new(j,i) >= eps_H) then
         mask(j,i) = 0_i1b   ! grounded ice
      else
         mask(j,i) = 1_i1b   ! ice-free land
      end if

#if (MARGIN==2)

   else   ! (mask(j,i) == 2_i1b); sea

      if (H_new(j,i) >= eps_H) then

         if ( H_new(j,i) < (rhosw_rho_ratio*H_sea_new(j,i)) ) then

#if (MARINE_ICE_FORMATION==1)
            mask(j,i) = 2_i1b   ! floating ice cut off -> sea
#elif (MARINE_ICE_FORMATION==2)
            mask(j,i) = 0_i1b   ! "underwater ice"
#endif
         else
            mask(j,i) = 0_i1b   ! grounded ice
         end if

#if (MARINE_ICE_CALVING==2 \
      || MARINE_ICE_CALVING==4 \
      || MARINE_ICE_CALVING==6)
         if (zl0(j,i) < z_mar) mask(j,i) = 2_i1b   ! sea
#elif (MARINE_ICE_CALVING==3 \
        || MARINE_ICE_CALVING==5 \
        || MARINE_ICE_CALVING==7)
         if (zl_new(j,i) < z_mar) mask(j,i) = 2_i1b   ! sea
#endif

      else
         mask(j,i) = 2_i1b   ! sea
      end if

#endif

   end if

end do
end do

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2)
if (time >= time_target_topo_final) mask = mask_target
#endif

!-------- Correction of zs_new and zb_new for ice-free land and sea --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 1_i1b) then   ! ice-free land

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 2_i1b) then   ! sea

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 3_i1b) then   ! floating ice

      errormsg = ' >>> calc_thk_mask_update_aux2:' &
               //           end_of_line &
               //'          mask(j,i)==3 not allowed for' &
               //           end_of_line &
               //'          hybrid SIA/SStA dynamics of ice sheets' &
               //           end_of_line &
               //'          without ice shelves!'
      call error(errormsg)

   end if

end do
end do

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice(z_sl)

!-------- Computation of further quantities --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 0_i1b) then   ! grounded ice

      if (n_cts(j,i) == 1_i1b) then
         if (H(j,i) > 0.0_dp) then
            H_inv        = 1.0_dp/H(j,i)
            H_c_new(j,i) = H_c(j,i) * H_new(j,i)*H_inv
            H_t_new(j,i) = H_t(j,i) * H_new(j,i)*H_inv
         else
            H_c_new(j,i) = 0.99_dp * H_new(j,i)   ! this case should not occur,
            H_t_new(j,i) = 0.01_dp * H_new(j,i)   ! just for safety
         end if
         zm_new(j,i) = zb_new(j,i)+H_t_new(j,i)
      else
         H_c_new(j,i) = H_new(j,i)
         H_t_new(j,i) = 0.0_dp
         zm_new(j,i)  = zb_new(j,i)
      end if

   else   ! mask(j,i) == 1_i1b or 2_i1b, ice-free land or sea

      H_c_new(j,i) = 0.0_dp
      H_t_new(j,i) = 0.0_dp
      H_new(j,i)   = 0.0_dp
      zs_new(j,i)  = zb_new(j,i)
      zm_new(j,i)  = zb_new(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux2

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for coupled SIA/SSA or
!! SIA/SStA/SSA dynamics of ice sheets with ice shelves.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux3(time, dtime, dtime_inv, z_sl, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dtime_inv, z_sl, z_mar

integer(i4b) :: i, j, kc, kt
real(dp)     :: rhosw_rho_ratio, rho_rhosw_ratio
real(dp)     :: H_inv
logical      :: flag_calving_event

real(dp), dimension(0:JMAX,0:IMAX) :: H_sea_new, H_balance

!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW

!-------- Update of the mask --------

zs_new = zb_new + H_new   ! ice surface topography
H_sea_new = z_sl - zl_new    ! sea depth
H_balance = H_sea_new*rhosw_rho_ratio

#if (THK_EVOL>=1)

do i=1, IMAX-1
do j=1, JMAX-1

   ! grounding_line migration check

   if ( ( mask(j,i) <= 1_i1b &
         .and.    (mask(j,i+1)>1_i1b.or.mask(j,i-1)>1_i1b &
               .or.mask(j+1,i)>1_i1b.or.mask(j-1,i)>1_i1b) ) &
        .or. &
        ( mask(j,i)>=2_i1b &
         .and.    (mask(j,i+1)==0_i1b.or.mask(j,i-1)==0_i1b  &
               .or.mask(j+1,i)==0_i1b.or.mask(j-1,i)==0_i1b) ) ) then

      if (H_new(j,i) >= eps_H) then

         if (H_new(j,i)<H_balance(j,i).and.zl_new(j,i)<z_sl) then
            mask(j,i)    = 3_i1b
            zb_new(j,i)   = z_sl-rho_rhosw_ratio*H_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         else
            mask(j,i)    = 0_i1b
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dzl_dtau(j,i)
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         end if

      else   ! if (H_new(j,i) <= eps_H) then

         if (zl_new(j,i)>z_sl) then
            mask(j,i)    = 1_i1b
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)
         else
            mask(j,i)    = 2_i1b
            zb_new(j,i)   = z_sl
            dzb_dtau(j,i) = 0.0_dp
            zs_new(j,i)   = z_sl
         end if

      end if

   else if (mask(j,i) <= 1_i1b) then   ! grounded ice or ice-free land

      if (H_new(j,i) >= eps_H) then   ! can change mask 0 or 1
         mask(j,i)    = 0_i1b
         zb_new(j,i)   = zl_new(j,i)
         dzb_dtau(j,i) = dzl_dtau(j,i)
         zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
      else
         mask(j,i)    = 1_i1b
         zb_new(j,i)   = zl_new(j,i)
         dzb_dtau(j,i) = dzl_dtau(j,i)
         zs_new(j,i)   = zl_new(j,i)
      end if

   else   ! if (mask(j,i)==2_i1b.or.mask(j,i)==3_i1b) then

      if (H_new(j,i) > eps_H) then

         if (H_new(j,i)<H_balance(j,i).and.zl_new(j,i)<z_sl) then
            mask(j,i)    = 3_i1b
            zb_new(j,i)   = z_sl-rho_rhosw_ratio*H_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         else
            mask(j,i)    = 0_i1b
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dzl_dtau(j,i)
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         end if

      else   ! if (H_new(j,i) <= eps_H) then

         if (zl_new(j,i)>z_sl) then
            mask(j,i)    = 1_i1b
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)
         else
            mask(j,i)    = 2_i1b
            zb_new(j,i)   = z_sl
            dzb_dtau(j,i) = 0.0_dp
            zs_new(j,i)   = z_sl
         end if

      end if

   end if

end do
end do

#if (ICE_SHELF_CALVING==2)

#if !defined(ALLOW_OPENAD) /* Normal */

do

   flag_calving_front_1 = .false.
   flag_calving_event   = .false.

   do i=1, IMAX-1
   do j=1, JMAX-1

      if ( (mask(j,i)==3_i1b) &   ! floating ice
           .and. &
             (    (mask(j,i+1)==2_i1b)   &   ! with
              .or.(mask(j,i-1)==2_i1b)   &   ! one
              .or.(mask(j+1,i)==2_i1b)   &   ! neighbouring
              .or.(mask(j-1,i)==2_i1b) ) &   ! sea point
         ) &
         flag_calving_front_1(j,i) = .true.   ! preliminary detection
                                              ! of the calving front

      if ( flag_calving_front_1(j,i).and.(H_new(j,i) < H_CALV) ) then
         flag_calving_event = .true.  ! calving event,
         mask(j,i)         = 2_i1b   ! floating ice point changes to sea point
      end if

   end do
   end do

   if (.not.flag_calving_event) exit

end do

#else /* OpenAD */

   flag_calving_event   = .true.

do while (flag_calving_event)

   flag_calving_front_1 = .false.
   flag_calving_event   = .false.

   do i=1, IMAX-1
   do j=1, JMAX-1

      if ( (mask(j,i)==3_i1b) &   ! floating ice
           .and. &
             (    (mask(j,i+1)==2_i1b)   &   ! with
              .or.(mask(j,i-1)==2_i1b)   &   ! one
              .or.(mask(j+1,i)==2_i1b)   &   ! neighbouring
              .or.(mask(j-1,i)==2_i1b) ) &   ! sea point
         ) &
         flag_calving_front_1(j,i) = .true.       ! preliminary detection
                                                  ! of the calving front

      if ( (flag_calving_front_1(j,i)).and.(H_new(j,i) < H_CALV) ) then
         flag_calving_event = .true.  ! calving event,
         mask(j,i)         = 2_i1b   ! floating ice point changes to sea point
      end if

   end do
   end do

end do

#endif /* Normal vs. OpenAD */

#elif (ICE_SHELF_CALVING==3)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i)==3_i1b) mask(j,i) = 2_i1b
      ! float-kill: all floating ice points changed to sea points

end do
end do

#elif (ICE_SHELF_CALVING==9)

#if (defined(MISMIPP))

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==3_i1b).and.(xi(i) > Lx)) then
      mask(j,i) = 2_i1b   ! floating ice point changes to sea point
   end if

end do
end do

#else

errormsg = ' >>> calc_thk_mask_update_aux3:' &
         //           end_of_line &
         //'          Option ICE_SHELF_CALVING==9' &
         //           end_of_line &
         //'          only defined for MISMIP+!'
call error(errormsg)

#endif

#endif

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2)
if (time >= time_target_topo_final) mask = mask_target
#endif

!-------- Correction of zs_new and zb_new
!         for ice-free land, sea and floating ice --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 1_i1b) then   ! ice-free land

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 2_i1b) then   ! sea

      zs_new(j,i) = z_sl          ! both zs_new(j,i) and zb_new(j,i)
      zb_new(j,i) = z_sl          ! set to the current sea level stand
      H_new(j,i)  = 0.0_dp

   else if (mask(j,i) == 3_i1b) then   ! floating ice

      H_new(j,i) = zs_new(j,i)-zb_new(j,i)   ! ice thickness

      zb_new(j,i) = z_sl - rho_rhosw_ratio*H_new(j,i)   ! floating condition
      zs_new(j,i) = zb_new(j,i) + H_new(j,i)

   end if

end do
end do

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice(z_sl)

!-------- Computation of further quantities --------

do i=0, IMAX
do j=0, JMAX

   if ( (mask(j,i) == 0_i1b).or.(mask(j,i) == 3_i1b) ) then
                                        ! grounded or floating ice

      if (n_cts(j,i) == 1_i1b) then
         if (H(j,i) > 0.0_dp) then
            H_inv        = 1.0_dp/H(j,i)
            H_c_new(j,i) = H_c(j,i) * H_new(j,i)*H_inv
            H_t_new(j,i) = H_t(j,i) * H_new(j,i)*H_inv
         else
            H_c_new(j,i) = 0.99_dp * H_new(j,i)   ! this case should not occur,
            H_t_new(j,i) = 0.01_dp * H_new(j,i)   ! just for safety
         end if
         zm_new(j,i) = zb_new(j,i)+H_t_new(j,i)
      else
         H_c_new(j,i) = H_new(j,i)
         H_t_new(j,i) = 0.0_dp
         zm_new(j,i)  = zb_new(j,i)
      end if

   else   ! mask(j,i) == 1_i1b or 2_i1b, ice-free land or sea

      H_c_new(j,i) = 0.0_dp
      H_t_new(j,i) = 0.0_dp
      H_new(j,i)   = 0.0_dp
      zs_new(j,i)  = zb_new(j,i)
      zm_new(j,i)  = zb_new(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux3

!-------------------------------------------------------------------------------
!> Limit thickness of isolated ice points.
!<------------------------------------------------------------------------------
subroutine limit_thickness_isolated_ice(z_sl)

implicit none

real(dp), intent(in) :: z_sl

integer(i4b) :: i, j
real(dp)     :: rho_rhosw_ratio
real(dp)     :: H_isol_max
logical      :: flag_kill_isolated_ice

rho_rhosw_ratio = RHO/RHO_SW

#if (defined(H_ISOL_MAX))
  H_isol_max = H_ISOL_MAX
#else
  H_isol_max = 1.0e+03_dp   ! default value of 1000 m
#endif

if (H_isol_max < eps_H) then
   H_isol_max = 0.0_dp
   flag_kill_isolated_ice = .true.
else
   flag_kill_isolated_ice = .false.
end if

do i=1, IMAX-1
do j=1, JMAX-1

   if (mask(j,i) == 0_i1b) then   ! grounded ice

      if ( ((mask(j,i+1) == 1_i1b).or.(mask(j,i+1) == 2_i1b)) &
           .and. ((mask(j,i-1) == 1_i1b).or.(mask(j,i-1) == 2_i1b)) &
           .and. ((mask(j+1,i) == 1_i1b).or.(mask(j+1,i) == 2_i1b)) &
           .and. ((mask(j-1,i) == 1_i1b).or.(mask(j-1,i) == 2_i1b)) &
         ) then   ! all nearest neighbours ice-free

         H_new(j,i)  = min(H_new(j,i), H_isol_max)
         zs_new(j,i) = zb_new(j,i)+H_new(j,i)

         if (flag_kill_isolated_ice) then
            if (zb_new(j,i) >= z_sl) then
               mask(j,i) = 1_i1b
            else
               mask(j,i) = 2_i1b
            end if
         end if

      end if

   else if (mask(j,i) == 3_i1b) then   ! floating ice

      if ( ((mask(j,i+1) == 1_i1b).or.(mask(j,i+1) == 2_i1b)) &
           .and. ((mask(j,i-1) == 1_i1b).or.(mask(j,i-1) == 2_i1b)) &
           .and. ((mask(j+1,i) == 1_i1b).or.(mask(j+1,i) == 2_i1b)) &
           .and. ((mask(j-1,i) == 1_i1b).or.(mask(j-1,i) == 2_i1b)) &
         ) then   ! all nearest neighbours ice-free

         H_new(j,i)  = min(H_new(j,i), H_isol_max)
         zb_new(j,i) = z_sl-rho_rhosw_ratio*H_new(j,i)
         zs_new(j,i) = zb_new(j,i)+H_new(j,i)

         if (flag_kill_isolated_ice) then
            mask(j,i) = 2_i1b
         end if

      end if

   end if

end do
end do

end subroutine limit_thickness_isolated_ice

!-------------------------------------------------------------------------------
!> Enforce connectivity of the ocean.
!<------------------------------------------------------------------------------
subroutine ocean_connect()

implicit none

integer(i4b)                           :: i, j
integer(i1b), dimension(0:JMAX,0:IMAX) :: mask_connect
integer(i1b), dimension(0:JMAX,0:IMAX) :: mask_connect_save, mask_connect_diff
logical                                :: flag_change

!-------- Determine connected area allowed to be ocean --------

mask_connect = 0_i1b

mask_connect(0:1         , :          ) = 1_i1b
mask_connect(JMAX-1:JMAX , :          ) = 1_i1b
mask_connect(:           , 0:1        ) = 1_i1b
mask_connect(:           , IMAX-1:IMAX) = 1_i1b

flag_change = .true.

do while (flag_change)

   mask_connect_save = mask_connect

   do i=1, IMAX-1
   do j=1, JMAX-1

      if (mask_connect_save(j,i) == 1_i1b) then
         if (mask(j  ,i+1) >= 2_i1b) mask_connect(j  ,i+1) = 1_i1b
         if (mask(j  ,i-1) >= 2_i1b) mask_connect(j  ,i-1) = 1_i1b
         if (mask(j+1,i  ) >= 2_i1b) mask_connect(j+1,i  ) = 1_i1b
         if (mask(j-1,i  ) >= 2_i1b) mask_connect(j-1,i  ) = 1_i1b
         if (mask(j+1,i+1) >= 2_i1b) mask_connect(j+1,i+1) = 1_i1b
         if (mask(j+1,i-1) >= 2_i1b) mask_connect(j+1,i-1) = 1_i1b
         if (mask(j-1,i+1) >= 2_i1b) mask_connect(j-1,i+1) = 1_i1b
         if (mask(j-1,i-1) >= 2_i1b) mask_connect(j-1,i-1) = 1_i1b
      end if

   end do
   end do

#if !defined(ALLOW_OPENAD) /* Normal */

   mask_connect_diff = abs(mask_connect-mask_connect_save)

   if (maxval(mask_connect_diff) > 0_i1b) then
      flag_change = .true.
   else
      flag_change = .false.
   end if

#else /* OpenAD */

   flag_change = .false.

   do i=0, IMAX
   do j=0, JMAX
      mask_connect_diff(j,i) = mask_connect(j,i)-mask_connect_save(j,i)
      if (mask_connect_diff(j,i) /= 0_i1b) flag_change = .true.
   end do
   end do

#endif /* Normal vs. OpenAD */

end do

!-------- Reset disconnected "ocean islands" to ice-free land --------

#if !defined(ALLOW_OPENAD) /* Normal */

where ((mask == 2_i1b).and.(mask_connect == 0_i1b)) mask = 1_i1b

#else /* OpenAD */

do i=0, IMAX
do j=0, JMAX
   if ((mask(j,i) == 2_i1b).and.(mask_connect(j,i) == 0_i1b)) &
      mask(j,i) = 1_i1b
end do
end do

#endif /* Normal vs. OpenAD */

end subroutine ocean_connect

!-------------------------------------------------------------------------------
!> Determination of the several components of the mass balance:
!! Runoff, calving, basal melt. 
!<------------------------------------------------------------------------------
subroutine account_mb_source(dtime, z_mar)

  implicit none

  real(dp), intent(in) :: dtime
  real(dp), intent(in) :: z_mar

  integer(i4b) :: i, j
  real(dp)     :: dtime_inv

  dtime_inv = 1.0_dp/dtime

!-------- Compute mass balance components --------

  mb_source_apl      = 0.0_dp
  as_perp_apl        = 0.0_dp
  accum_apl          = 0.0_dp
  runoff_apl         = 0.0_dp
  Q_b_apl            = 0.0_dp
  calving_apl        = 0.0_dp
  mask_ablation_type = 0_i1b

  do i=0, IMAX
  do j=0, JMAX

     if (flag_inner_point(j,i)) then   ! inner point

        if (H_new(j,i) >= eps_H) then   ! glaciated point

           mb_source_apl(j,i) = (H_new(j,i)-H_new_flow(j,i))*dtime_inv
                                ! applied MB, based on ice thickness change

           ! Store melt quantities here for later accounting
           accum_apl(j,i)   = accum(j,i)
           runoff_apl(j,i)  = runoff(j,i)
           Q_b_apl(j,i)     = Q_b_tot(j,i)
           calving_apl(j,i) = calving(j,i)
           as_perp_apl(j,i) = accum_apl(j,i) - runoff_apl(j,i)

           ! Store melting mask value. This became an ice point,
           ! should be either grounded or floating ice.
           if (mask(j,i)==0_i1b) then
              mask_ablation_type(j,i) = 1_i1b ! grounded ice
           else if (mask(j,i)==3_i1b) then
              mask_ablation_type(j,i) = 3_i1b ! floating ice
           else
              mask_ablation_type(j,i) = 9_i1b ! misaccounted (if appears)
           end if

        else if (H_new(j,i) < eps_H .and. H_new_flow(j,i) >= eps_H) then
                                           ! ice disappeared

           mb_source_apl(j,i) = (H_new(j,i)-H_new_flow(j,i))*dtime_inv
                                ! applied MB, based on ice thickness change

           accum_apl(j,i) = accum(j,i)

           if (mask(j,i)==2_i1b) then
              ! all mass that flows into the ocean counted as
              ! (large-scale) calving
              calving_apl(j,i) = -(mb_source_apl(j,i)-accum_apl(j,i))

           else if (runoff(j,i) > 0.0_dp .or. calving(j,i) > 0.0_dp &
                                         .or. Q_b_tot(j,i) > 0.0_dp) then
              ! land or shelf where ice disappered,
              ! three melting contingents shared proportionally
              runoff_apl(j,i)  = -(mb_source_apl(j,i)-accum_apl(j,i)) &
                                       *runoff(j,i) &
                                       /(runoff(j,i)+calving(j,i) &
                                                    +Q_b_tot(j,i))
              calving_apl(j,i) = -(mb_source_apl(j,i)-accum_apl(j,i)) &
                                       *calving(j,i) &
                                       /(runoff(j,i)+calving(j,i) &
                                                    +Q_b_tot(j,i))
              Q_b_apl(j,i)     = -(mb_source_apl(j,i)-accum_apl(j,i)) &
                                       *Q_b_tot(j,i) &
                                       /(runoff(j,i)+calving(j,i) &
                                                    +Q_b_tot(j,i))

           !!! else runoff(j,i)=0.0_dp .and. calving(j,i) = 0.0_dp &
           !!!                         .and. Q_b_tot(j,i) = 0.0_dp

           end if

           as_perp_apl(j,i) = accum_apl(j,i) - runoff_apl(j,i)

           ! Store melting mask value.
           ! This grid point is ice free now, therefore it
           ! can be either ocean or land (using actual mask value).

           if (mask(j,i) == 2_i1b) then   ! ocean point (hidden melt)
              mask_ablation_type(j,i) = -2_i1b
           else                            ! land (hidden melt)
              mask_ablation_type(j,i) = -1_i1b
           end if

        end if

     end if

  end do
  end do

end subroutine account_mb_source

!-------------------------------------------------------------------------------

end module calc_thk_m
!
