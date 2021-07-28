
MODULE SICO_MATHS_M_DIFF
  USE SICO_TYPES_M
  IMPLICIT NONE
!-------------------------------------------------------------------------------
  PUBLIC 

CONTAINS

subroutine transpose_csr(a_value, a_index, a_diag_index, a_ptr, &
                    nnz, nmax, &
                    b_value, b_index, b_diag_index, b_ptr)

implicit none

integer(i4b),                     intent(in)  :: nnz, nmax
integer(i4b), dimension(nmax+1),  intent(in)  :: a_ptr
integer(i4b), dimension(nnz),     intent(in)  :: a_index
integer(i4b), dimension(nmax),    intent(in)  :: a_diag_index
real(dp),     dimension(nnz),     intent(in)  :: a_value
integer(i4b), dimension(nmax+1),  intent(out) :: b_ptr
integer(i4b), dimension(nnz),     intent(out) :: b_index
integer(i4b), dimension(nmax),    intent(out) :: b_diag_index
real(dp),     dimension(nnz),     intent(out) :: b_value

integer(i4b) :: nr, k
integer(i4b), dimension(nmax)                 :: b_ptr_plus

   b_ptr = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_ptr(a_index(k)+1) = b_ptr(a_index(k)+1) + 1
      end do
   end do
   b_ptr(1) = 1
   do nr=1, nmax
     b_ptr(nr+1) = b_ptr(nr) + b_ptr(nr+1)
   end do
   b_ptr_plus = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_index(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = nr
         if(a_index(k).eq.nr)then
           b_diag_index(nr) = b_ptr(a_index(k)) + b_ptr_plus(a_index(k))
         end if
         b_value(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = a_value(k)
         b_ptr_plus(a_index(k)) = b_ptr_plus(a_index(k)) + 1
      end do
   end do
end subroutine transpose_csr

subroutine transpose_csr_no_diagonal(a_value, a_index, a_ptr, &
                    nnz, nmax, &
                    b_value, b_index, b_ptr)

implicit none

integer(i4b),                     intent(in)  :: nnz, nmax
integer(i4b), dimension(nmax+1),  intent(in)  :: a_ptr
integer(i4b), dimension(nnz),     intent(in)  :: a_index
real(dp),     dimension(nnz),     intent(in)  :: a_value
integer(i4b), dimension(nmax+1),  intent(out) :: b_ptr
integer(i4b), dimension(nnz),     intent(out) :: b_index
real(dp),     dimension(nnz),     intent(out) :: b_value

integer(i4b) :: nr, k
integer(i4b), dimension(nmax)                 :: b_ptr_plus

   b_ptr = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_ptr(a_index(k)+1) = b_ptr(a_index(k)+1) + 1
      end do
   end do
   b_ptr(1) = 1
   do nr=1, nmax
     b_ptr(nr+1) = b_ptr(nr) + b_ptr(nr+1)
   end do
   b_ptr_plus = 0
   do nr=1, nmax
      do k=a_ptr(nr), a_ptr(nr+1)-1
         b_index(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = nr
         b_value(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = a_value(k)
         b_ptr_plus(a_index(k)) = b_ptr_plus(a_index(k)) + 1
      end do
   end do
end subroutine transpose_csr_no_diagonal


!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
subroutine sor_sprs_b(lgs_a_value, lgs_a_valueb, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                    lgs_b_value, lgs_b_valueb, &
                    nnz, nmax, n_sprs, omega, eps_sor, lgs_x_value, lgs_x_valueb, ierr)

implicit none

integer(i4b),                     intent(in) :: nnz, nmax, n_sprs
real(dp),                         intent(in) :: omega, eps_sor
integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
integer(i4b), dimension(nnz),     intent(in) :: lgs_a_index
integer(i4b), dimension(nmax),    intent(in) :: lgs_a_diag_index
real(dp),     dimension(nnz),     intent(in) :: lgs_a_value
real(dp),     dimension(nmax),    intent(in) :: lgs_b_value

integer(i4b),                    intent(out) :: ierr
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

real(dp),     dimension(nnz),  intent(inout) :: lgs_a_valueb
real(dp),     dimension(nmax), intent(inout) :: lgs_b_valueb
real(dp),     dimension(nmax), intent(inout) :: lgs_x_valueb

integer(i4b), dimension(nmax+1)              :: lgs_aT_ptr
integer(i4b), dimension(nnz)                 :: lgs_aT_index
integer(i4b), dimension(nmax)                :: lgs_aT_diag_index
real(dp),     dimension(nnz)                 :: lgs_aT_value
real(dp),     dimension(nmax)                :: incrbb
integer(i4b) :: iter
integer(i4b) :: nr, k
real(dp), allocatable, dimension(:) :: lgs_x_value_prev
real(dp)     :: b_nr
logical      :: flag_convergence
logical      :: inif, infor

  incrbb = 0.0
  call transpose_csr(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                    nnz, nmax, &
                    lgs_aT_value, lgs_aT_index, lgs_aT_diag_index, lgs_aT_ptr)

  call sor_sprs(lgs_aT_value, lgs_aT_index, lgs_aT_diag_index, lgs_aT_ptr, &
                    lgs_x_valueb, &
                    nnz, nmax, n_sprs, omega, eps_sor, incrbb, ierr)

  DO nr=1,nmax
    lgs_b_valueb(nr) = lgs_b_valueb(nr) + incrbb(nr)
  ENDDO

  call sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                    lgs_b_value, &
                    nnz, nmax, n_sprs, omega, eps_sor, lgs_x_value, ierr)
  do nr=1, nmax
    inif = .false. 
    infor = .false. 
    if(nr .eq. lgs_a_diag_index(nr)) then
      inif = .true.
      lgs_a_valueb(lgs_a_diag_index(nr)) = lgs_a_valueb(lgs_a_diag_index(nr)) -&
             lgs_x_value(lgs_a_index(nr)) * incrbb(nr)
    end if
    do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      infor = .true.
      lgs_a_valueb(k) = lgs_a_valueb(k) - lgs_x_value(lgs_a_index(k)) * incrbb(nr)
    end do
  end do
  lgs_x_valueb = 0.0
end subroutine sor_sprs_b

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
  SUBROUTINE SOR_SPRS(lgs_a_value, lgs_a_index, lgs_a_diag_index, &
&   lgs_a_ptr, lgs_b_value, nnz, nmax, n_sprs, omega, eps_sor, &
&   lgs_x_value, ierr)
  implicit none

  integer(i4b),                     intent(in) :: n_sprs
  integer(i4b),                     intent(in) :: nnz, nmax
  real(dp),                         intent(in) :: omega, eps_sor
  integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),     intent(in) :: lgs_a_index
  integer(i4b), dimension(nmax),    intent(in) :: lgs_a_diag_index
  real(dp),     dimension(nnz),     intent(in) :: lgs_a_value
  real(dp),     dimension(nmax),    intent(in) :: lgs_b_value

  integer(i4b),                    intent(out) :: ierr
  real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

  integer(i4b) :: iter
  integer(i4b) :: iter_max
  integer(i4b) :: nr, k
  real(dp),           dimension(nmax) :: lgs_x_value_prev
  real(dp)     :: temp1, temp2
  logical      :: isnanflag1, isnanflag2, isnanflag3
  real(dp)     :: b_nr
  logical      :: flag_convergence

#if (ITER_MAX_SOR > 0)
  iter_max = ITER_MAX_SOR
#else
  iter_max = 1000   ! default value
#endif


  iter_loop : do iter=1, iter_max

     lgs_x_value_prev = lgs_x_value

     do nr=1, nmax

        b_nr = 0.0_dp 

        do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
           b_nr = b_nr + lgs_a_value(k)*lgs_x_value(lgs_a_index(k))
        end do

        lgs_x_value(nr) = lgs_x_value(nr) &
                          -omega*(b_nr-lgs_b_value(nr)) &
                                /lgs_a_value(lgs_a_diag_index(nr))

     end do

     flag_convergence = .true.
     do nr=1, nmax
        if (abs(lgs_x_value(nr)-lgs_x_value_prev(nr)) > eps_sor) then
           flag_convergence = .false.
           exit
        end if
     end do

     if (flag_convergence) then
        write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
        ierr = 0   ! convergence criterion fulfilled
        return
     end if

  end do iter_loop

  write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
  ierr = -1   ! convergence criterion not fulfilled

  END SUBROUTINE SOR_SPRS

!  Differentiation of tri_sle in reverse (adjoint) mode:
!   gradient     of useful results: x b
!   with respect to varying inputs: x a0 a1 a2 b
!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!! @param[in]  a0       a0(j) is element A_(j,j-1) of Matrix A
!! @param[in]  a1       a1(j) is element A_(j,j)   of Matrix A
!! @param[in]  a2       a2(j) is element A_(j,j+1) of Matrix A
!! @param[in]  b        inhomogeneity vector
!! @param[in]  nrows    size of matrix A (indices run from 0 (!!!) to nrows)
!! @param[out] x        Solution vector.
!<------------------------------------------------------------------------------
  SUBROUTINE TRI_SLE_B(a0, a0b, a1, a1b, a2, a2b, x, xb, b, bb, nrows)
    IMPLICIT NONE
    INTEGER(i4b), INTENT(IN) :: nrows
    REAL(dp), DIMENSION(0:nrows), INTENT(IN) :: a0, a2
    REAL(dp), DIMENSION(0:nrows) :: a0b, a2b
    REAL(dp), DIMENSION(0:nrows), INTENT(INOUT) :: a1, b
    REAL(dp), DIMENSION(0:nrows), INTENT(INOUT) :: a1b, bb
    REAL(dp), DIMENSION(0:nrows) :: x
    REAL(dp), DIMENSION(0:nrows) :: xb
    REAL(dp), DIMENSION(0:nrows) :: help_x
    REAL(dp), DIMENSION(0:nrows) :: help_xb
    INTEGER(i4b) :: n
    REAL(dp) :: tempb
    REAL(dp) :: tempb0
!--------  Generate an upper triangular matrix
!                      ('obere Dreiecksmatrix') --------
    DO n=1,nrows
      CALL PUSHREAL8(a1(n))
      a1(n) = a1(n) - a0(n)/a1(n-1)*a2(n-1)
    END DO
    DO n=1,nrows
      CALL PUSHREAL8(b(n))
      b(n) = b(n) - a0(n)/a1(n-1)*b(n-1)
    END DO
    help_x(0) = b(nrows)/a1(nrows)
    DO n=1,nrows
      CALL PUSHREAL8(help_x(n))
      help_x(n) = b(nrows-n)/a1(nrows-n) - a2(nrows-n)/a1(nrows-n)*&
&       help_x(n-1)
    END DO
    help_xb = 0.0_8
    DO n=nrows,0,-1
      help_xb(nrows-n) = help_xb(nrows-n) + xb(n)
      xb(n) = 0.0_8
    END DO
    DO n=nrows,1,-1
      CALL POPREAL8(help_x(n))
      tempb = help_xb(n)/a1(nrows-n)
      tempb0 = -(help_xb(n)/a1(nrows-n))
      help_xb(n) = 0.0_8
      a2b(nrows-n) = a2b(nrows-n) + help_x(n-1)*tempb0
      help_xb(n-1) = help_xb(n-1) + a2(nrows-n)*tempb0
      a1b(nrows-n) = a1b(nrows-n) - a2(nrows-n)*help_x(n-1)*tempb0/a1(&
&       nrows-n) - b(nrows-n)*tempb/a1(nrows-n)
      bb(nrows-n) = bb(nrows-n) + tempb
    END DO
    tempb = help_xb(0)/a1(nrows)
    bb(nrows) = bb(nrows) + tempb
    a1b(nrows) = a1b(nrows) - b(nrows)*tempb/a1(nrows)
    DO n=nrows,1,-1
      CALL POPREAL8(b(n))
      tempb = -(bb(n)/a1(n-1))
      a0b(n) = a0b(n) + b(n-1)*tempb
      bb(n-1) = bb(n-1) + a0(n)*tempb
      a1b(n-1) = a1b(n-1) - a0(n)*b(n-1)*tempb/a1(n-1)
    END DO
    DO n=nrows,1,-1
      CALL POPREAL8(a1(n))
      tempb = -(a1b(n)/a1(n-1))
      a0b(n) = a0b(n) + a2(n-1)*tempb
      a2b(n-1) = a2b(n-1) + a0(n)*tempb
      a1b(n-1) = a1b(n-1) - a0(n)*a2(n-1)*tempb/a1(n-1)
    END DO
  END SUBROUTINE TRI_SLE_B

!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!! @param[in]  a0       a0(j) is element A_(j,j-1) of Matrix A
!! @param[in]  a1       a1(j) is element A_(j,j)   of Matrix A
!! @param[in]  a2       a2(j) is element A_(j,j+1) of Matrix A
!! @param[in]  b        inhomogeneity vector
!! @param[in]  nrows    size of matrix A (indices run from 0 (!!!) to nrows)
!! @param[out] x        Solution vector.
!<------------------------------------------------------------------------------
  SUBROUTINE TRI_SLE(a0, a1, a2, x, b, nrows)
    IMPLICIT NONE
!       (The trick with the help_x was introduced in order to avoid
!        the negative step in the original, blanked-out loop.)
!deallocate(help_x)
!  WARNING: Subroutine does not check for elements of the main
!           diagonal becoming zero. In this case it crashes even
!           though the system may be solvable. Otherwise ok.
    INTEGER(i4b), INTENT(IN) :: nrows
    REAL(dp), DIMENSION(0:nrows), INTENT(IN) :: a0, a2
    REAL(dp), DIMENSION(0:nrows), INTENT(INOUT) :: a1, b
    REAL(dp), DIMENSION(0:nrows), INTENT(OUT) :: x
    REAL(dp), DIMENSION(0:nrows) :: help_x
    INTEGER(i4b) :: n
!--------  Generate an upper triangular matrix
!                      ('obere Dreiecksmatrix') --------
    DO n=1,nrows
      a1(n) = a1(n) - a0(n)/a1(n-1)*a2(n-1)
    END DO
    DO n=1,nrows
      b(n) = b(n) - a0(n)/a1(n-1)*b(n-1)
    END DO
! a0(n)  = 0.0_dp , not needed in the following, therefore
!                   not set
!-------- Iterative solution of the new system --------
! x(nrows) = b(nrows)/a1(nrows)
! do n=nrows-1, 0, -1
!    x(n) = (b(n)-a2(n)*x(n+1))/a1(n)
! end do
!allocate(help_x(0:nrows))
    help_x(0) = b(nrows)/a1(nrows)
    DO n=1,nrows
      help_x(n) = b(nrows-n)/a1(nrows-n) - a2(nrows-n)/a1(nrows-n)*&
&       help_x(n-1)
    END DO
    DO n=0,nrows
      x(n) = help_x(nrows-n)
    END DO
  END SUBROUTINE TRI_SLE

!-------------------------------------------------------------------------------
!> Bilinear interpolation.
!<------------------------------------------------------------------------------
  FUNCTION BILININT(x1, x2, y1, y2, z11, z12, z21, z22, x, y)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x1, x2, y1, y2, z11, z12, z21, z22, x, y
    REAL(dp) :: t, u
    REAL(dp) :: bilinint
    REAL(dp), PARAMETER :: i=1.0_dp
    t = (x-x1)/(x2-x1)
    u = (y-y1)/(y2-y1)
    bilinint = (i-t)*(i-u)*z11 + (i-t)*u*z12 + t*(i-u)*z21 + t*u*z22
  END FUNCTION BILININT

!  Differentiation of my_erfc in reverse (adjoint) mode:
!   gradient     of useful results: retval
!   with respect to varying inputs: x
!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!<------------------------------------------------------------------------------
  SUBROUTINE MY_ERFC_B(x, xb, retval, retvalb)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x
    REAL(dp) :: xb
    REAL(dp) :: retval
    REAL(dp) :: retvalb
    REAL(dp) :: t, z
    REAL(dp) :: tb, zb
    INTRINSIC ABS
    INTRINSIC EXP
    REAL(dp) :: temp
    REAL(dp) :: temp0
    REAL(dp) :: temp1
    REAL(dp) :: temp2
    REAL(dp) :: temp3
    REAL(dp) :: temp4
    REAL(dp) :: tempb
    REAL(dp) :: tempb0
    REAL(dp) :: tempb1
    REAL(dp) :: tempb2
    INTEGER :: branch
    IF (x .GE. 0.) THEN
      z = x
      CALL PUSHCONTROL1B(0)
    ELSE
      z = -x
      CALL PUSHCONTROL1B(1)
    END IF
    t = 1.0_dp/(1.0_dp+0.5_dp*z)
    IF (x .LT. 0.0_dp) retvalb = -retvalb
    temp = t*(0.17087277_dp*t-0.82215223_dp) + 1.48851587_dp
    temp0 = t*temp - 1.13520398_dp
    temp1 = t*(t*temp0+0.27886807_dp) - 0.18628806_dp
    temp2 = t*(t*temp1+0.09678418_dp) + 0.37409196_dp
    temp3 = t*temp2 + 1.00002368_dp
    temp4 = t*temp3 - z*z - 1.26551223_dp
    tempb = EXP(temp4)*t*retvalb
    tempb0 = t**2*tempb
    tempb1 = t**2*tempb0
    tempb2 = t**2*tempb1
    tb = EXP(temp4)*retvalb + (temp3+temp2*t)*tempb + (t*temp1+temp1*t+&
&     0.09678418_dp)*tempb0 + (t*temp0+temp0*t+0.27886807_dp)*tempb1 + (&
&     temp+(0.17087277_dp*t-0.82215223_dp)*t+0.17087277_dp*t**2)*tempb2
    zb = -(2*z*tempb) - 0.5_dp*tb/(0.5_dp*z+1.0_dp)**2
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xb = zb
    ELSE
      xb = -zb
    END IF
  END SUBROUTINE MY_ERFC_B

!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!<------------------------------------------------------------------------------
  SUBROUTINE MY_ERFC(x, retval)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(OUT) :: retval
    REAL(dp) :: t, z
    INTRINSIC ABS
    INTRINSIC EXP
    IF (x .GE. 0.) THEN
      z = x
    ELSE
      z = -x
    END IF
    t = 1.0_dp/(1.0_dp+0.5_dp*z)
    retval = t*EXP(-(z*z)-1.26551223_dp+t*(1.00002368_dp+t*(&
&     0.37409196_dp+t*(0.09678418_dp+t*(-0.18628806_dp+t*(0.27886807_dp+&
&     t*(-1.13520398_dp+t*(1.48851587_dp+t*(-0.82215223_dp+t*&
&     0.17087277_dp)))))))))
    IF (x .LT. 0.0_dp) retval = 2.0_dp - retval
  END SUBROUTINE MY_ERFC

#if  defined(BUILD_LIS) && (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
!#if  (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
!-------------------------------------------------------------------------------
!> OpenAD needs a template to help differentiate through the LIS solver when  
!! it is used. This is substituted in for adjoint modes in Antarctica with 
!! ice shelves.
!<------------------------------------------------------------------------------
#include "lisf.h"
  SUBROUTINE SICO_LIS_SOLVER(nmax, nnz, lgs_a_ptr, lgs_a_index, &
&   lgs_a_value, lgs_b_value, lgs_x_value)
implicit none

integer(i4b)                                 :: ierr
integer(i4b)                                 :: iter
integer(i4b)                                 :: nc, nr
integer(i4b),                     intent(in) :: nmax
integer(i4b),                     intent(in) :: nnz
integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
integer(i4b), dimension(nnz),  intent(in)    :: lgs_a_index

LIS_MATRIX                                   :: lgs_a
LIS_VECTOR                                   :: lgs_b
LIS_VECTOR                                   :: lgs_x
LIS_SOLVER                                   :: solver

real(dp),     dimension(nnz),  intent(in)    :: lgs_a_value
real(dp),     dimension(nmax),    intent(in) :: lgs_b_value
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

character(len=256)                           :: ch_solver_set_option

!  ------ Settings for Lis
call lis_init_f(ierr)
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

#if (defined(LIS_OPTS))
    ch_solver_set_option = trim(LIS_OPTS)
#else
    ch_solver_set_option = '-i bicgsafe -p jacobi '// &
                           '-maxiter 2000 -tol 1.0e-18 -initx_zeros false'
#endif

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
call lis_finalize_f(ierr)
  END SUBROUTINE SICO_LIS_SOLVER

subroutine sico_lis_solver_b(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value, lgs_a_valueb, lgs_b_value,& 
                           lgs_b_valueb, lgs_x_value, lgs_x_valueb)

implicit none

integer(i4b)                                 :: nc, nr
integer(i4b),                     intent(in) :: nmax
integer(i4b),                     intent(in) :: nnz
integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
integer(i4b), dimension(nnz),  intent(in) :: lgs_a_index

real(dp),     dimension(nnz),  intent(in) :: lgs_a_value
real(dp),     dimension(nmax),    intent(in) :: lgs_b_value
real(dp),     dimension(nmax), intent(inout) :: lgs_x_value

real(dp),     dimension(nnz),  intent(inout) :: lgs_a_valueb
real(dp),     dimension(nmax), intent(inout) :: lgs_b_valueb
real(dp),     dimension(nmax), intent(inout) :: lgs_x_valueb

integer(i4b), dimension(nmax+1)              :: lgs_aT_ptr
integer(i4b), dimension(nnz)              :: lgs_aT_index
real(dp),     dimension(nnz)              :: lgs_aT_value
real(dp),     dimension(nmax)                :: incrbb

lgs_aT_value = 0.
lgs_aT_index = 0
lgs_aT_value = 0.
  call transpose_csr_no_diagonal(lgs_a_value, lgs_a_index, lgs_a_ptr, &
                    nnz, nmax, &
                    lgs_aT_value, lgs_aT_index, lgs_aT_ptr)

incrbb = 0.
print *, "incrbb before: ", SUM(ABS(incrbb))
  call sico_lis_solver(nmax, nnz, &
                           lgs_aT_ptr, lgs_aT_index, &
                           lgs_aT_value, lgs_x_valueb, incrbb)
print *, "incrbb after: ", SUM(ABS(incrbb))
  DO nr=1,nmax
    lgs_b_valueb(nr) = lgs_b_valueb(nr) + incrbb(nr)
  ENDDO

  call sico_lis_solver(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value, lgs_b_value, lgs_x_value)

  do nr=1, nmax
    do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      lgs_a_valueb(nc) = lgs_a_valueb(nc) - lgs_x_value(lgs_a_index(nc)) * incrbb(nr)
    end do
  end do
  lgs_x_valueb = 0.0
end subroutine sico_lis_solver_b

#endif
END MODULE SICO_MATHS_M_DIFF
!

