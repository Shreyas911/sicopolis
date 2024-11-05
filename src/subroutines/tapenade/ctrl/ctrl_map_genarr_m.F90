module ctrl_map_genarr_m

  use sico_types_m
  use sico_variables_m
  use error_m

  implicit none

#ifdef DO_CTRL_GENARR2D
  public :: ctrl_map_ini_genarr2d, ctrl_map_genarr2d
#endif
#ifdef DO_CTRL_GENARR3D
  public :: ctrl_map_ini_genarr3d, ctrl_map_genarr3d
#endif

#if (defined(DO_CTRL_GENARR2D) || defined(DO_CTRL_GENARR3D))
contains
#endif

#ifdef DO_CTRL_GENARR2D
  subroutine ctrl_map_ini_genarr2d

    implicit none

    integer(i4b)        :: ctrl_index
    integer(i4b)        :: igen_c_slide_init, igen_H, igen_q_geo
    integer(i4b)        :: igen_gamma_s, igen_s_stat
    integer(i4b)        :: igen_beta1, igen_beta2, igen_c_dis_da
    integer(i4b)        :: igen_Pmax, igen_mu, igen_delta_tda_const

#ifdef XX_GENARR2D_VARS_ARR
    xx_genarr2d_vars            = XX_GENARR2D_VARS_ARR
#endif
#ifdef XX_GENARR2D_PREPROC_ARR
    xx_genarr2d_preproc         = XX_GENARR2D_PREPROC_ARR
#endif
#ifdef XX_GENARR2D_LOG10INITVAL_ARR
    xx_genarr2d_log10initval    = XX_GENARR2D_LOG10INITVAL_ARR
#endif

    igen_c_slide_init    = 0
    igen_H               = 0
    igen_q_geo           = 0
    igen_gamma_s         = 0
    igen_s_stat          = 0
    igen_beta1           = 0
    igen_beta2           = 0
    igen_Pmax            = 0
    igen_mu              = 0
    igen_delta_tda_const = 0
    igen_c_dis_da        = 0

#ifdef XX_GENARR2D_VARS_ARR
    do ctrl_index = 1, NUM_CTRL_GENARR2D
      
      if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_c_slide_init') then
        igen_c_slide_init = ctrl_index
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_q_geo') then
        igen_q_geo = ctrl_index
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_H') then
        igen_H = ctrl_index
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_gamma_s') then
        igen_gamma_s = ctrl_index
#if (ACCSURFACE != 2 && ACCSURFACE != 3)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'gamma_s as a control param is only compatible with ACCSURFACE == 2 or ACCSURFACE == 3 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_s_stat') then
        igen_s_stat = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'s_stat as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 '&
        //'or (ACCSURFACE <= 5 && SOLID_PRECIP == 3)!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_beta1') then
        igen_beta1 = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'beta1 as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_beta2') then
        igen_beta2 = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'beta2 as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_Pmax') then
        igen_Pmax = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'Pmax as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_mu') then
        igen_mu = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'mu as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_delta_tda_const') then
        igen_delta_tda_const = ctrl_index
#if (TSURFACE > 4)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'delta_tda_const as a control param is only compatible with TSURFACE <= 4 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_c_dis_da') then
        igen_c_dis_da = ctrl_index
#if (!defined(GRL) || DISC <= 0)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'c_dis_da as a control param is only compatible with GRL domain and DISC > 0!'
#endif
      else
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
          //"This control variable is not in the genctrl2d setup yet!"
        call error(errormsg)
      end if
  
    end do

    if (igen_c_slide_init .GT. 0) then
      call ctrl_map_genarr2d(c_slide_init, igen_c_slide_init)
    end if
    if (igen_q_geo .GT. 0) then
      call ctrl_map_genarr2d(q_geo, igen_q_geo)
    end if
    if (igen_H .GT. 0) then
      call ctrl_map_genarr2d(H, igen_H)
    end if
#if (ACCSURFACE==2 || ACCSURFACE==3)
    if (igen_gamma_s .GT. 0) then
      call ctrl_map_genarr2d(gamma_s, igen_gamma_s)
    end if
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2 || (ACCSURFACE<=5 && SOLID_PRECIP==3))
    if (igen_s_stat .GT. 0) then
      call ctrl_map_genarr2d(s_stat, igen_s_stat)
    end if
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2)
    if (igen_beta1 .GT. 0) then
      call ctrl_map_genarr2d(beta1, igen_beta1)
    end if
    if (igen_beta2 .GT. 0) then
      call ctrl_map_genarr2d(beta2, igen_beta2)
    end if
    if (igen_Pmax .GT. 0) then
      call ctrl_map_genarr2d(Pmax, igen_Pmax)
    end if
    if (igen_mu .GT. 0) then
      call ctrl_map_genarr2d(mu, igen_mu)
    end if
    if (igen_delta_tda_const .GT. 0) then
      call ctrl_map_genarr2d(delta_tda_const, igen_delta_tda_const)
    end if
    if (igen_c_dis_da .GT. 0) then
      call ctrl_map_genarr2d(c_dis_da, igen_c_dis_da)
    end if
#endif
#endif
    
  end subroutine ctrl_map_ini_genarr2d

  subroutine ctrl_map_genarr2d(fld, iarr)

    implicit none  

    real(dp), dimension(0:JMAX,0:IMAX)          :: fld
    integer(i4b)                                :: iarr, k2
    logical                                     :: dopreprocs, dolog10ctrl
    real(dp)                                    :: log10initval, ln10
#ifdef XX_GENARR2D_PREPROC_ARR
    character(128), dimension(NUMCTRLPROCARR2D) :: preprocs
#else
! Dummy definition
    character(128), dimension(1) :: preprocs
#endif

    dopreprocs = .FALSE.
    ln10 = log(10.0)
    dolog10ctrl = .FALSE.

#ifdef XX_GENARR2D_PREPROC_ARR
    dopreprocs = .TRUE.
#endif

    if (dopreprocs) then

      read (unit=xx_genarr2d_preproc(iarr),fmt=*) preprocs

      do k2 = 1, NUMCTRLPROCARR2D
#ifdef XX_GENARR2D_LOG10INITVAL_ARR
        if (preprocs(k2) .EQ. 'log10ctrl') then
          dolog10ctrl = .TRUE.
          log10initval = xx_genarr2d_log10initval(iarr)
        end if
#endif
      end do

#ifdef XX_GENARR2D_LOG10INITVAL_ARR
      if (dolog10ctrl) then  
        xx_genarr2d(iarr,:,:) = xx_genarr2d(iarr,:,:) + log10initval
        xx_genarr2d(iarr,:,:) = EXP(ln10 * xx_genarr2d(iarr,:,:)) 
        fld = xx_genarr2d(iarr,:,:)
      else
        fld = fld + xx_genarr2d(iarr,:,:)
      endif
#endif

    else

      fld = fld + xx_genarr2d(iarr,:,:)

    endif

  end subroutine ctrl_map_genarr2d
#endif

#ifdef DO_CTRL_GENARR3D
  subroutine ctrl_map_ini_genarr3d

    implicit none

    integer(i4b)        :: ctrl_index
    integer(i4b)        :: igen_vx_c, igen_vy_c, igen_vz_c
    integer(i4b)        :: igen_temp_c, igen_age_c

#ifdef XX_GENARR3D_VARS_ARR
    xx_genarr3d_vars            = XX_GENARR3D_VARS_ARR
#endif
#ifdef XX_GENARR3D_PREPROC_ARR
    xx_genarr3d_preproc         = XX_GENARR3D_PREPROC_ARR
#endif
#ifdef XX_GENARR3D_LOG10INITVAL_ARR
    xx_genarr3d_log10initval    = XX_GENARR3D_LOG10INITVAL_ARR
#endif

    igen_vx_c = 0
    igen_vy_c = 0
    igen_vz_c = 0
    igen_temp_c = 0
    igen_age_c = 0

#ifdef XX_GENARR3D_VARS_ARR
    do ctrl_index = 1, NUM_CTRL_GENARR3D
      if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_vx_c') then
        igen_vx_c = ctrl_index
#if (!(ANF_DAT==3 && RESTART==1))
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //'vx_c as a control param is only compatible with ' &
          //'ANF_DAT == 3 and RESTART == 1 for now !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_vy_c') then
        igen_vy_c = ctrl_index
#if (!(ANF_DAT==3 && RESTART==1))
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //'vy_c as a control param is only compatible with ' &
          //'ANF_DAT == 3 and RESTART == 1 for now !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_vz_c') then
        igen_vz_c = ctrl_index
#if (!(ANF_DAT==3 && RESTART==1))
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //'vz_c as a control param is only compatible with ' &
          //'ANF_DAT == 3 and RESTART == 1 for now !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_temp_c') then
        igen_temp_c = ctrl_index
#if (ANF_DAT==1 && TEMP_INIT==5)
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //'Initial ice temperature as a control param is incompatible with ' &
          //'ANF_DAT == 1 and TEMP_INIT == 5 for now !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_age_c') then
        igen_age_c = ctrl_index       
      else
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //"This control variable is not in the genctrl3d setup yet!"
        call error(errormsg)
      end if
  
    end do

    if (igen_vx_c .GT. 0) then
      call ctrl_map_genarr3d(vx_c, igen_vx_c)
    end if
    if (igen_vy_c .GT. 0) then
      call ctrl_map_genarr3d(vy_c, igen_vy_c)
    end if
    if (igen_vz_c .GT. 0) then
      call ctrl_map_genarr3d(vz_c, igen_vz_c)
    end if
    if (igen_temp_c .GT. 0) then
      call ctrl_map_genarr3d(temp_c, igen_temp_c)
    end if
    if (igen_age_c .GT. 0) then
      call ctrl_map_genarr3d(age_c, igen_age_c)
    end if
#endif
    
  end subroutine ctrl_map_ini_genarr3d

  subroutine ctrl_map_genarr3d(fld, iarr)

    implicit none  

    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)  :: fld
    integer(i4b)                                :: iarr, k3
    logical                                     :: dopreprocs, dolog10ctrl
    real(dp)                                    :: log10initval, ln10
#ifdef XX_GENARR3D_PREPROC_ARR
    character(128), dimension(NUMCTRLPROCARR3D) :: preprocs
#else
! Dummy definition
    character(128), dimension(1) :: preprocs
#endif

    dopreprocs = .FALSE.
    ln10 = log(10.0)
    dolog10ctrl = .FALSE.

#ifdef XX_GENARR3D_PREPROC_ARR
    dopreprocs = .TRUE.
#endif

    if (dopreprocs) then

      read (unit=xx_genarr3d_preproc(iarr),fmt=*) preprocs

      do k3 = 1, NUMCTRLPROCARR3D
#ifdef XX_GENARR3D_LOG10INITVAL_ARR
        if (preprocs(k3) .EQ. 'log10ctrl') then
          dolog10ctrl = .TRUE.
          log10initval = xx_genarr3d_log10initval(iarr)
        end if
#endif
      end do

#ifdef XX_GENARR3D_LOG10INITVAL_ARR
      if (dolog10ctrl) then  
        xx_genarr3d(iarr,:,:,:) = xx_genarr3d(iarr,:,:,:) + log10initval
        xx_genarr3d(iarr,:,:,:) = EXP(ln10 * xx_genarr3d(iarr,:,:,:))
        fld = xx_genarr3d(iarr,:,:,:)
      else
        fld = fld + xx_genarr3d(iarr,:,:,:)
      endif
#endif

    else

      fld = fld + xx_genarr3d(iarr,:,:,:)

    end if       

  end subroutine ctrl_map_genarr3d
#endif

end module ctrl_map_genarr_m 
