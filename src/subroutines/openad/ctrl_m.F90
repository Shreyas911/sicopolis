!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!   Module  :  c t r l
  
!   Purpose :  Declarations of control variables for adjointing

!   Copyright 2017-2022 Liz Curry-Logan, Sri Hari Krishna Narayanan
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
module ctrl_m

  use sico_types_m  
  use sico_variables_m
  use sico_vars_m

  implicit none

  public :: ctrl_init
  public :: cost_independent_init, cost_dependent_init
  public :: cost_final 
#if (defined(ALLOW_GRDCHK) || defined(ALLOW_OPENAD))
  public :: myceiling, myfloor 
#endif

contains
 
!------------------------------------------------------------------------------- 
!> Initialiation of control variable.
!! Recognized by OpenAD with the prefix xxVar,
!! where Var is the variable of choice (normally
!! something in sico_variables_m
!<------------------------------------------------------------------------------
  subroutine ctrl_init()
  
  implicit none

  integer(i4b) :: i, j, k 

  ! CHOICE OF CONTROL VARIABLE IS MADE HERE:
  double precision, dimension(0:JMAX,0:IMAX) ::         xxH_c 
  double precision, dimension(0:JMAX,0:IMAX) ::         xxc_drag
  double precision, dimension(0:JMAX,0:IMAX) ::         xxc_slide
  double precision, dimension(0:JMAX,0:IMAX) ::         xxtemp_ma_present
  double precision, dimension(0:KCMAX,0:JMAX,0:IMAX) :: xxtemp_c
  double precision, dimension(0:JMAX,0:IMAX) ::         xxq_geo
  double precision, dimension(0:JMAX,0:IMAX) ::         xxvis_int_g
  double precision, dimension(0:JMAX,0:IMAX) ::         xxdzs_dxi_g
  double precision, dimension(0:JMAX,0:IMAX) ::         xxdzs_deta_g
  double precision, dimension(0:JMAX,0:IMAX) ::         xxcalving
  double precision, dimension(0:JMAX,0:IMAX) ::         xxQ_bm
  double precision, dimension(0:JMAX,0:IMAX) ::         xxQ_tld
  double precision, dimension(0:JMAX,0:IMAX) ::         xxacc_fact
  double precision, dimension(0:JMAX,0:IMAX,12) ::      xxprecip_present
  double precision, dimension(0:KCMAX,0:JMAX,0:IMAX) :: xxsigma_c
  double precision, dimension(0:KCMAX,0:JMAX,0:IMAX) :: xxvx_c

  ! 3D controls:
  do k=0, KCMAX
    do j=0, JMAX
      do i=0, IMAX
        xxtemp_c(k,j,i) = 0.0d0
        xxvx_c(k,j,i) = 0.0d0
        xxsigma_c(k,j,i) = 0.0d0
      end do
    end do
  end do
 
  ! 2D controls:
  do i=0, IMAX
    do j=0, JMAX
  
      ! IC: ice geometry
      xxH_c(j,i) = 0.0d0
  
      ! BC: sliding parameter 
      xxc_drag(j,i) = 0.0d0
      xxc_slide(j,i) = 0.0d0
  
      ! BC: temperature-related quantities:
      xxq_geo(j,i) = 0.0d0
      xxtemp_ma_present(j,i) = 0.0d0
      xxvis_int_g(j,i) = 0.0d0
      xxdzs_dxi_g(j,i) = 0.0d0
      xxdzs_deta_g(j,i) = 0.0d0
      xxcalving(j,i) = 0.0d0
      xxQ_bm(j,i) = 0.0d0
      xxQ_tld(j,i) = 0.0d0

      ! BC: SMB-related quantities 
      xxprecip_present(j,i,1) = 0.0d0
      xxacc_fact(j,i) = 0.0d0
 
    end do
  end do
  
  end subroutine ctrl_init

!-------------------------------------------------------------------------------
!> The independent variable is perturbed by the xxVar
!! at all places during the adjoint mode.
!<------------------------------------------------------------------------------
  subroutine cost_independent_init()
  
  implicit none
  
  integer(i4b) :: i, j, k 
  double precision, dimension(0:JMAX,0:IMAX) ::         xxH_c 
  double precision, dimension(0:JMAX,0:IMAX) ::         xxtemp_ma_present
  double precision, dimension(0:KCMAX,0:JMAX,0:IMAX) :: xxtemp_c
  double precision, dimension(0:KCMAX,0:JMAX,0:IMAX) :: xxsigma_c
  double precision, dimension(0:KCMAX,0:JMAX,0:IMAX) :: xxvx_c
  double precision, dimension(0:JMAX,0:IMAX,12) ::      xxprecip_present 
  double precision, dimension(0:JMAX,0:IMAX) ::         xxacc_fact
  double precision, dimension(0:JMAX,0:IMAX) ::         xxq_geo
  double precision, dimension(0:JMAX,0:IMAX) ::         xxc_drag
  double precision, dimension(0:JMAX,0:IMAX) ::         xxc_slide
  double precision, dimension(0:JMAX,0:IMAX) ::         xxvis_int_g
  double precision, dimension(0:JMAX,0:IMAX) ::         xxdzs_dxi_g
  double precision, dimension(0:JMAX,0:IMAX) ::         xxdzs_deta_g
  double precision, dimension(0:JMAX,0:IMAX) ::         xxcalving
  double precision, dimension(0:JMAX,0:IMAX) ::         xxQ_bm
  double precision, dimension(0:JMAX,0:IMAX) ::         xxQ_tld
  
  !-------- Initialize independent ctrl variables --------

  ! 3D controls:
  do k=0, KCMAX
    do j=0, JMAX
      do i=0, IMAX
        temp_c(k,j,i) = temp_c(k,j,i) + xxtemp_c(k,j,i)
        sigma_c(k,j,i) = sigma_c(k,j,i) + xxsigma_c(k,j,i)
        vx_c(k,j,i) = vx_c(k,j,i) + xxvx_c(k,j,i)
      end do
    end do
  end do

  ! 2D controls:
  do j=0,JMAX
      do i=0,IMAX
 
        ! IC: ice geometry
        H_c(j,i) = H_c(j,i) + xxH_c(j,i)
  
        ! Basal processes: 
        c_drag(j,i) = c_drag(j,i) + xxc_drag(j,i)
        c_slide(j,i) = c_slide(j,i) + xxc_slide(j,i)
 
        ! BC: temperature-related quantities 
        q_geo(j,i) = q_geo(j,i) + xxq_geo(j,i)
        temp_ma_present(j,i) = temp_ma_present(j,i) + xxtemp_ma_present(j,i)
        vis_int_g(j,i) = vis_int_g(j,i) + xxvis_int_g(j,i)
        dzs_dxi_g(j,i) = dzs_dxi_g(j,i) + xxdzs_dxi_g(j,i)
        dzs_deta_g(j,i) = dzs_deta_g(j,i) + xxdzs_deta_g(j,i)
        calving(j,i) = calving(j,i) + xxcalving(j,i)
        Q_bm(j,i) = Q_bm(j,i) + xxQ_bm(j,i)
        Q_tld(j,i) = Q_tld(j,i) + xxQ_tld(j,i)
  
        ! BC: SMB-related quantities 
        precip_present(j,i,1) = precip_present(j,i,1) + xxprecip_present(j,i,1)
        acc_fact(j,i) = acc_fact(j,i) + xxacc_fact(j,i)
 
      end do
  end do
  
  end subroutine cost_independent_init
  
!-------------------------------------------------------------------------------
!> The dependent variable for the cost routine.
!! This must be a scalar variable, although option for
!! a summation or weighted average of multiple cost 
!! targets is possible.
!<------------------------------------------------------------------------------  
  subroutine cost_dependent_init()
  
  implicit none
  
  !------- Final cost value:
  fc = 0.0_dp
  
  !------- Can be sum of multiple cost test values (later):
  objf_test = 0.0_dp
  mult_test = 1.0_dp
  
  end subroutine cost_dependent_init

!-------------------------------------------------------------------------------
!> This is the final cost calculation.
!! The cost function structure is defined here. Currently 
!! is a "observed age" - modeled age summed over the entire 
!! domain. The "observed age" is a fake, generated age
!! field performed by the 125 ka run in headers.
!! Other cost functions (e.g., total cold ice volume,
!! commented out below) are certainly possible, and 
!! recommended! 
!<------------------------------------------------------------------------------
  subroutine cost_final(runname)
  
  implicit none
  
  integer(i4b) :: i, j, k, kc, kt, ios, KDATA 
  character(len=100), intent(out) :: runname
  
  !-------- Calculate the difference between the modeled and 'observed' ages:
#ifdef AGE_COST

#if (CALCMOD!=1)
  KDATA = KCMAX
#else 
  KDATA = KCMAX + KTMAX
print *, '>>> error: CALCMOD == 1 but final cost not properly working for '
print *, '           AGE_COST simulations'
#endif
 
  do k=0, KDATA 
    do i=0, IMAX
      do j=0, JMAX

        ! only counting points that are real in the data: 
        if (  age_data(k,j,i).lt. 0.0       &
         .and.age_unc(k,j,i) .lt. 0.0    ) then

          objf_test = objf_test + &
                      sqrt(((age_data(k,j,i) - age_c(k,j,i)))**2 &
                      / ((age_unc(k,j,i))**2))

        end if
      end do
    end do
  end do

#else

    do i=0, IMAX
      do j=0, JMAX
        !--- Other cost functions:
        objf_test = objf_test + (H_c(j,i) + H_t(j,i))*cell_area(j,i)
      end do
    end do

#endif

  fc = mult_test * (objf_test)
  
  !-------- Print to screen just in case something gets
  !         crazy with the file outputting:
  print *, 'Final cost, fc = ', fc
  print *, trim(OUT_PATH)
  
  !-------- Write final cost to a file:
  open(unit=97, iostat=ios, &
#ifndef ALLOW_OPENAD
       file=trim(OUT_PATH)//'/'//trim(runname)//'_COST.dat', &
#else
       file='AD_COST', &
#endif
       status='new')
  write(unit=97, fmt='(f26.6)') fc
  write(unit=97, fmt='(2a)') 'Final cost, fc = ', runname
  close(unit=97)
  
  end subroutine cost_final

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_OPENAD))
  subroutine myceiling(num, output)
  
  implicit none
  
  integer(i4b)             :: inum
  real(dp), intent(in) :: num
  integer(i4b), intent(out) :: output
  
  inum = int(num);
  
  inum = int(num);
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    output = inum;
  else if (num>0) then
    output = inum + 1;
  else if (num<0) then
    output = inum;
  end if
  end subroutine myceiling
  
  subroutine myfloor(num, output)
  
  implicit none
  
  integer(i4b)             :: inum
  real(dp), intent(in) :: num
  integer(i4b), intent(out) :: output
  
  inum = int(num);
  
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    output = inum;
  else if (num>0) then
    output = inum;
  else if (num<0) then
    output = inum - 1;
  end if
  
  end subroutine myfloor
#endif

end module ctrl_m
