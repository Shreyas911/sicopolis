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
#ifdef DO_CTRL_GENARR3DR
  public :: ctrl_map_ini_genarr3dr, ctrl_map_genarr3dr
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
    integer(i4b)        :: igen_RHO_A, igen_time_lag_asth, igen_flex_rig_lith
    integer(i4b)        :: igen_p_weert, igen_q_weert
    integer(i4b)        :: igen_enh_fact_da_dummy2d_scalar, igen_enh_intg_da_dummy2d_scalar
    integer(i4b)        :: igen_n_glen_da_dummy2d_scalar
    integer(i4b)        :: igen_zs, igen_zl, igen_zl0, igen_zb

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
    igen_RHO_A           = 0
    igen_time_lag_asth   = 0
    igen_flex_rig_lith   = 0
    igen_p_weert         = 0
    igen_q_weert         = 0
    igen_enh_fact_da_dummy2d_scalar = 0
    igen_enh_intg_da_dummy2d_scalar = 0
    igen_n_glen_da_dummy2d_scalar   = 0
    igen_zs              = 0
    igen_zl              = 0
    igen_zl0             = 0
    igen_zb              = 0

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
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_RHO_A') then
        igen_RHO_A = ctrl_index
#if (!defined(PARAM_RHO_A))
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'RHO_A as a control param is only compatible with PARAM_RHO_A being defined!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_time_lag_asth') then
        igen_time_lag_asth = ctrl_index
#if (REBOUND != 1 && REBOUND != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'time_lag_asth as a control param is only compatible with REBOUND == 1 or REBOUND == 2!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_flex_rig_lith') then
        igen_flex_rig_lith = ctrl_index
#if (REBOUND != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'flex_rig_lith as a control param is only compatible with REBOUND == 2!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_p_weert') then
        igen_p_weert = ctrl_index
#if (!defined(P_WEERT))
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'p_weert as a control param is only compatible with P_WEERT being defined!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_q_weert') then
        igen_q_weert = ctrl_index
#if (!defined(Q_WEERT))
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'q_weert as a control param is only compatible with Q_WEERT being defined!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_enh_fact_da_dummy2d_scalar') then
        igen_enh_fact_da_dummy2d_scalar = ctrl_index
#if (ENHMOD != 1 && ENHMOD != 2 && ENHMOD != 3)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'enh_fact_da_dummy2d_scalar as a control param is only compatible with ENHMOD == 1, 2, 3!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_enh_intg_da_dummy2d_scalar') then
        igen_enh_intg_da_dummy2d_scalar = ctrl_index
#if (ENHMOD != 2 && ENHMOD != 3)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'enh_intg_da_dummy2d_scalar as a control param is only compatible with ENHMOD == 2, 3!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_n_glen_da_dummy2d_scalar') then
        igen_n_glen_da_dummy2d_scalar = ctrl_index
#if (FLOW_LAW != 1)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'n_glen_da_dummy2d_scalar as a control param is only compatible with FLOW_LAW == 1!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_zs') then
        igen_zs = ctrl_index
#if (ANF_DAT==2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'zs as a control param is only compatible with ANF_DAT != 2!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_zl') then
        igen_zl = ctrl_index
#if (ANF_DAT==2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'zl as a control param is only compatible with ANF_DAT != 2!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_zl0') then
        igen_zl0 = ctrl_index
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_zb') then
        igen_zb = ctrl_index
#if (ANF_DAT==2 || (ANF_DAT == 1 && !defined(ZB_PRESENT_FILE)))
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'zb as a control param is only compatible with ANF_DAT == 1 and defined(ZB_PRESENT_FILE) or ANF_DAT == 3!'
        call error(errormsg)
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
#endif
#if (defined(GRL) && DISC > 0)
    if (igen_c_dis_da .GT. 0) then
      call ctrl_map_genarr2d(c_dis_da, igen_c_dis_da)
    end if
#endif
#if (defined(PARAM_RHO_A))
    if (igen_RHO_A .GT. 0) then
      call ctrl_map_genarr2d(RHO_A, igen_RHO_A)
    end if
#endif
#if (REBOUND==1 || REBOUND==2)
    if (igen_time_lag_asth .GT. 0) then
      call ctrl_map_genarr2d(time_lag_asth, igen_time_lag_asth)
    end if
#endif
#if (REBOUND==2)
    if (igen_flex_rig_lith .GT. 0) then
      call ctrl_map_genarr2d(flex_rig_lith, igen_flex_rig_lith)
    end if
#endif
#if (defined(P_WEERT))
    if (igen_p_weert .GT. 0) then
      call ctrl_map_genarr2d(p_weert, igen_p_weert)
    end if
#endif
#if (defined(Q_WEERT))
    if (igen_q_weert .GT. 0) then
      call ctrl_map_genarr2d(q_weert, igen_q_weert)
    end if
#endif
#if (ENHMOD==1 || ENHMOD==2 || ENHMOD==3)
    if (igen_enh_fact_da_dummy2d_scalar .GT. 0) then
      call ctrl_map_genarr2d(enh_fact_da_dummy2d_scalar, igen_enh_fact_da_dummy2d_scalar)
    end if
#endif
#if (ENHMOD==2 || ENHMOD==3)
    if (igen_enh_intg_da_dummy2d_scalar .GT. 0) then
      call ctrl_map_genarr2d(enh_intg_da_dummy2d_scalar, igen_enh_intg_da_dummy2d_scalar)
    end if
#endif
#if (FLOW_LAW==1)
    if (igen_n_glen_da_dummy2d_scalar .GT. 0) then
      call ctrl_map_genarr2d(n_glen_da_dummy2d_scalar, igen_n_glen_da_dummy2d_scalar)
    end if
#endif
#if (ANF_DAT!=2)
    if (igen_zs .GT. 0) then
      call ctrl_map_genarr2d(zs, igen_zs)
    end if
    if (igen_zl .GT. 0) then
      call ctrl_map_genarr2d(zl, igen_zl)
    end if
#endif
    if (igen_zl0 .GT. 0) then
      call ctrl_map_genarr2d(zl0, igen_zl0)
    end if
#if ((ANF_DAT==1 && defined(ZB_PRESENT_FILE)) || ANF_DAT==3)
    if (igen_zb .GT. 0) then
      call ctrl_map_genarr2d(zb, igen_zb)
    end if
#endif

#endif
    
  end subroutine ctrl_map_ini_genarr2d

  subroutine ctrl_map_genarr2d(fld, iarr)

    implicit none  

    real(dp), dimension(0:JMAX,0:IMAX)          :: fld
    integer(i4b)                                :: iarr, k2
    logical                                     :: dolog10ctrl
    real(dp)                                    :: log10initval, ln10
#ifdef XX_GENARR2D_PREPROC_ARR
    character(128), dimension(NUMCTRLPROCARR2D) :: preprocs
#endif

    ln10 = log(10.0)
    dolog10ctrl = .FALSE.

#ifdef XX_GENARR2D_PREPROC_ARR

    read (unit=xx_genarr2d_preproc(iarr),fmt=*) preprocs

    do k2 = 1, NUMCTRLPROCARR2D
      if (preprocs(k2) .EQ. 'log10ctrl') then
        dolog10ctrl = .TRUE.
#if (defined(XX_GENARR2D_LOG10INITVAL_ARR) && !defined(AD_INPUT_PATH))
        log10initval = xx_genarr2d_log10initval(iarr)
#endif
      end if
    end do

    if (dolog10ctrl) then
#if (defined(XX_GENARR2D_LOG10INITVAL_ARR) && !defined(AD_INPUT_PATH))
      xx_genarr2d(iarr,:,:) = xx_genarr2d(iarr,:,:) + log10initval
#endif
      xx_genarr2d(iarr,:,:) = EXP(ln10 * xx_genarr2d(iarr,:,:))
      fld = xx_genarr2d(iarr,:,:)
    else
      fld = fld + xx_genarr2d(iarr,:,:)
    endif

#else
    fld = fld + xx_genarr2d(iarr,:,:)
#endif

  end subroutine ctrl_map_genarr2d
#endif

#ifdef DO_CTRL_GENARR3D
  subroutine ctrl_map_ini_genarr3d

    implicit none

    integer(i4b)        :: ctrl_index
    integer(i4b)        :: igen_omega_c, igen_temp_c, igen_age_c

#ifdef XX_GENARR3D_VARS_ARR
    xx_genarr3d_vars            = XX_GENARR3D_VARS_ARR
#endif
#ifdef XX_GENARR3D_PREPROC_ARR
    xx_genarr3d_preproc         = XX_GENARR3D_PREPROC_ARR
#endif
#ifdef XX_GENARR3D_LOG10INITVAL_ARR
    xx_genarr3d_log10initval    = XX_GENARR3D_LOG10INITVAL_ARR
#endif

    igen_omega_c = 0
    igen_temp_c = 0
    igen_age_c = 0

#ifdef XX_GENARR3D_VARS_ARR
    do ctrl_index = 1, NUM_CTRL_GENARR3D
      if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_omega_c') then
        igen_omega_c = ctrl_index
#if !(ANF_DAT==3)
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //'omega_c as a control param is only compatible with ' &
          //'ANF_DAT == 3 for now !'
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

#if (ANF_DAT==3)
    if (igen_omega_c .GT. 0) then
      call ctrl_map_genarr3d(omega_c, igen_omega_c)
    end if
#endif
#if ((ANF_DAT!=1) || ((ANF_DAT==1) && (TEMP_INIT!=5)))
    if (igen_temp_c .GT. 0) then
      call ctrl_map_genarr3d(temp_c, igen_temp_c)
    end if
#endif
    if (igen_age_c .GT. 0) then
      call ctrl_map_genarr3d(age_c, igen_age_c)
    end if
#endif
    
  end subroutine ctrl_map_ini_genarr3d

  subroutine ctrl_map_genarr3d(fld, iarr)

    implicit none  

    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)  :: fld
    integer(i4b)                                :: iarr, k3
    logical                                     :: dolog10ctrl
    real(dp)                                    :: log10initval, ln10
#ifdef XX_GENARR3D_PREPROC_ARR
    character(128), dimension(NUMCTRLPROCARR3D) :: preprocs
#endif

    ln10 = log(10.0)
    dolog10ctrl = .FALSE.

#ifdef XX_GENARR3D_PREPROC_ARR

    read (unit=xx_genarr3d_preproc(iarr),fmt=*) preprocs

    do k3 = 1, NUMCTRLPROCARR3D
      if (preprocs(k3) .EQ. 'log10ctrl') then
        dolog10ctrl = .TRUE.
#if (defined(XX_GENARR3D_LOG10INITVAL_ARR) && !defined(AD_INPUT_PATH))
        log10initval = xx_genarr3d_log10initval(iarr)
#endif
      end if
    end do

    if (dolog10ctrl) then
#if (defined(XX_GENARR3D_LOG10INITVAL_ARR) && !defined(AD_INPUT_PATH))
      xx_genarr3d(iarr,:,:,:) = xx_genarr3d(iarr,:,:,:) + log10initval
#endif
      xx_genarr3d(iarr,:,:,:) = EXP(ln10 * xx_genarr3d(iarr,:,:,:))
      fld = xx_genarr3d(iarr,:,:,:)
    else
      fld = fld + xx_genarr3d(iarr,:,:,:)
    endif

#else
    fld = fld + xx_genarr3d(iarr,:,:,:)
#endif

  end subroutine ctrl_map_genarr3d
#endif

#ifdef DO_CTRL_GENARR3DR
  subroutine ctrl_map_ini_genarr3dr

    implicit none

    integer(i4b)        :: ctrl_index
    integer(i4b)        :: igen_temp_r

#ifdef XX_GENARR3DR_VARS_ARR
    xx_genarr3dr_vars            = XX_GENARR3DR_VARS_ARR
#endif
#ifdef XX_GENARR3DR_PREPROC_ARR
    xx_genarr3dr_preproc         = XX_GENARR3DR_PREPROC_ARR
#endif
#ifdef XX_GENARR3DR_LOG10INITVAL_ARR
    xx_genarr3dr_log10initval    = XX_GENARR3DR_LOG10INITVAL_ARR
#endif

    igen_temp_r = 0

#ifdef XX_GENARR3DR_VARS_ARR
    do ctrl_index = 1, NUM_CTRL_GENARR3DR
      if (trim(adjustl(xx_genarr3dr_vars(ctrl_index))) .EQ. 'xx_temp_r') then
        igen_temp_r = ctrl_index
#if !(ANF_DAT==3)
        errormsg = ' >>> ctrl_map_ini_genarr3dr: ' &
          //'temp_r as a control param is only compatible with ' &
          //'ANF_DAT == 3 for now !'
        call error(errormsg)
#endif
      else
        errormsg = ' >>> ctrl_map_ini_genarr3dr: ' &
          //"This control variable is not in the genctrl3dr setup yet!"
        call error(errormsg)
      end if

    end do

#if (ANF_DAT==3)
    if (igen_temp_r .GT. 0) then
      call ctrl_map_genarr3dr(temp_r, igen_temp_r)
    end if
#endif
#endif

  end subroutine ctrl_map_ini_genarr3dr

  subroutine ctrl_map_genarr3dr(fld, iarr)

    implicit none

    real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX)  :: fld
    integer(i4b)                                :: iarr, k3
    logical                                     :: dolog10ctrl
    real(dp)                                    :: log10initval, ln10
#ifdef XX_GENARR3DR_PREPROC_ARR
    character(128), dimension(NUMCTRLPROCARR3DR) :: preprocs
#endif

    ln10 = log(10.0)
    dolog10ctrl = .FALSE.

#ifdef XX_GENARR3DR_PREPROC_ARR

    read (unit=xx_genarr3dr_preproc(iarr),fmt=*) preprocs

    do k3 = 1, NUMCTRLPROCARR3DR
      if (preprocs(k3) .EQ. 'log10ctrl') then
        dolog10ctrl = .TRUE.
#if (defined(XX_GENARR3DR_LOG10INITVAL_ARR) && !defined(AD_INPUT_PATH))
        log10initval = xx_genarr3dr_log10initval(iarr)
#endif
      end if
    end do

    if (dolog10ctrl) then
#if (defined(XX_GENARR3DR_LOG10INITVAL_ARR) && !defined(AD_INPUT_PATH))
      xx_genarr3dr(iarr,:,:,:) = xx_genarr3dr(iarr,:,:,:) + log10initval
#endif
      xx_genarr3dr(iarr,:,:,:) = EXP(ln10 * xx_genarr3dr(iarr,:,:,:))
      fld = xx_genarr3dr(iarr,:,:,:)
    else
      fld = fld + xx_genarr3dr(iarr,:,:,:)
    endif

#else
    fld = fld + xx_genarr3dr(iarr,:,:,:)
#endif

  end subroutine ctrl_map_genarr3dr
#endif

end module ctrl_map_genarr_m 
