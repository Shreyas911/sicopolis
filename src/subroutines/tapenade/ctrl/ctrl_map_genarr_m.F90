#if defined(ALLOW_GENCTRL)
module ctrl_map_genarr_m

  use sico_types_m
  use sico_variables_m
  use error_m

  implicit none

  public :: ctrl_map_ini_genarr2d, ctrl_map_genarr2d
  public :: ctrl_map_ini_genarr3d, ctrl_map_genarr3d

contains

  subroutine ctrl_map_ini_genarr2d

    implicit none

    integer(i4b)        :: ctrl_index
    integer(i4b)        :: igen_c_slide_init, igen_H, igen_q_geo
    integer(i4b)        :: igen_gamma_s_arr, igen_s_stat_arr
    integer(i4b)        :: igen_beta1_arr_orig, igen_beta2_arr_orig
    integer(i4b)        :: igen_Pmax_arr, igen_mu_arr_orig

    xx_genarr2d_vars            = XX_GENARR2D_VARS_ARR
    xx_genarr2d_preproc         = XX_GENARR2D_PREPROC_ARR
    xx_genarr2d_bounds          = XX_GENARR2D_BOUNDS_ARR
    xx_genarr2d_log10initval    = XX_GENARR2D_LOG10INITVAL_ARR
    xx_genarr2d_weight          = XX_GENARR2D_WEIGHT_ARR

    igen_c_slide_init   = 0
    igen_H              = 0
    igen_q_geo          = 0
    igen_gamma_s_arr    = 0
    igen_s_stat_arr     = 0
    igen_beta1_arr_orig = 0
    igen_beta2_arr_orig = 0
    igen_Pmax_arr       = 0
    igen_mu_arr_orig    = 0

    do ctrl_index = 1, NUM_CTRL_GENARR2D
      
      if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_c_slide_init') then
        igen_c_slide_init = ctrl_index
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_q_geo') then
        igen_q_geo = ctrl_index
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_H') then
        igen_H = ctrl_index
#if (ANF_DAT != 1)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
          //'Initial surface topography as a control param is only compatible with ANF_DAT == 1 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_gamma_s_arr') then
        igen_gamma_s_arr = ctrl_index
#if (ACCSURFACE != 2 && ACCSURFACE != 3)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'gamma_s_arr as a control param is only compatible with ACCSURFACE == 2 or ACCSURFACE == 3 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_s_stat_arr') then
        igen_s_stat_arr = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2 && (ACCSURFACE > 5 || SOLID_PRECIP != 3))
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'s_stat_arr as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 '&
        //'or (ACCSURFACE <= 5 && SOLID_PRECIP == 3)!'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_beta1_arr_orig') then
        igen_beta1_arr_orig = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'beta1_arr_orig as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_beta2_arr_orig') then
        igen_beta2_arr_orig = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'beta2_arr_orig as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_Pmax_arr') then
        igen_Pmax_arr = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'Pmax_arr as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
        call error(errormsg)
#endif
      else if (trim(adjustl(xx_genarr2d_vars(ctrl_index))) .EQ. 'xx_mu_arr_orig') then
        igen_mu_arr_orig = ctrl_index
#if (ABLSURFACE != 1 && ABLSURFACE != 2)
        errormsg = ' >>> ctrl_map_ini_genarr2d: ' &
        //'mu_arr_orig as a control param is only compatible with ABLSURFACE == 1 or ABLSURFACE == 2 !'
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
#if (ANF_DAT == 1)
    if (igen_H .GT. 0) then
      call ctrl_map_genarr2d(H, igen_H)
    end if
#endif
#if (ACCSURFACE==2 || ACCSURFACE==3)
    if (igen_gamma_s_arr .GT. 0) then
      call ctrl_map_genarr2d(gamma_s_arr, igen_gamma_s_arr)
    end if
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2 || (ACCSURFACE<=5 && SOLID_PRECIP==3))
    if (igen_s_stat_arr .GT. 0) then
      call ctrl_map_genarr2d(s_stat_arr, igen_s_stat_arr)
    end if
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2)
    if (igen_beta1_arr_orig .GT. 0) then
      call ctrl_map_genarr2d(beta1_arr_orig, igen_beta1_arr_orig)
    end if
    if (igen_beta2_arr_orig .GT. 0) then
      call ctrl_map_genarr2d(beta2_arr_orig, igen_beta2_arr_orig)
    end if
    if (igen_Pmax_arr .GT. 0) then
      call ctrl_map_genarr2d(Pmax_arr, igen_Pmax_arr)
    end if
    if (igen_mu_arr_orig .GT. 0) then
      call ctrl_map_genarr2d(mu_arr_orig, igen_mu_arr_orig)
    end if
#endif
    
  end subroutine ctrl_map_ini_genarr2d

  subroutine ctrl_map_genarr2d(fld, iarr)

    implicit none  

    real(dp), dimension(0:JMAX,0:IMAX)          :: fld, wgenarr2d
    integer(i4b)                                :: iarr, j, i, k2
    integer(i4b)                                :: iLow, iHigh, jLow, jHigh
    logical                                     :: doread, dobounds
    logical                                     :: dosmooth, doscaling, dolog10ctrl
    real(dp)                                    :: log10initval, ln10
    character(128), dimension(NUMCTRLPROCARR2D)      :: preprocs

    ln10 = log(10.0)

    doread      = .FALSE.
    dosmooth    = .FALSE.
    doscaling   = .FALSE.
    dolog10ctrl = .FALSE.
    dobounds    = .FALSE.

    read (unit=xx_genarr2d_preproc(iarr),fmt=*) preprocs

    do k2 = 1, NUMCTRLPROCARR2D

      if (preprocs(k2) .EQ. 'smooth') then
        dosmooth = .TRUE.
        errormsg = ' >>> ctrl_map_genarr2d: ' &
          //'Weaver and Courtier like smoothing is not yet implemented!'
        call error(errormsg)
      end if
      if (preprocs(k2) .EQ. 'scaling') then
        doscaling = .TRUE.
        do i = 0, IMAX
          do j = 0, JMAX
            read(unit=xx_genarr2d_weight(iarr), fmt=*) wgenarr2d(j,i)
          end do
        end do
      else
        wgenarr2d = 1.0
      end if
      if (preprocs(k2) .EQ. 'read') then
        doread = .TRUE.
        errormsg = ' >>> ctrl_map_genarr2d: ' &
          //'Reading xx_gen from a file is yet to be implemented!'
        call error(errormsg)
      end if
      if (preprocs(k2) .EQ. 'log10ctrl') then
        dolog10ctrl = .TRUE.
        log10initval = xx_genarr2d_log10initval(iarr)
      end if
      if (preprocs(k2) .EQ. 'bounds') then
        dobounds = .TRUE.
      end if
    end do

    if (dobounds .AND. trim(adjustl(xx_genarr2d_bounds(iarr))) .NE. ' ') then
      read (unit=xx_genarr2d_bounds(iarr),fmt=*) iLow, iHigh, jLow, jHigh
    else
      iLow  = 0
      iHigh = IMAX
      jLow  = 0
      jHigh = JMAX
    end if

    do i = 0, IMAX
      do j = 0, JMAX
        if ((i .GE. iLow) .AND. (i .LE. iHigh) .AND. (j .GE. jLow) .AND. (j .LE. jHigh)) then
          xx_genarr2d_mask(iarr,j,i) = 1.0
        else
          xx_genarr2d_mask(iarr,j,i) = 0.0
        end if
      end do
    end do

    xx_genarr2d(iarr,:,:) = xx_genarr2d(iarr,:,:) / sqrt(wgenarr2d(:,:))

    if (dolog10ctrl) then  
      xx_genarr2d(iarr,:,:) = xx_genarr2d(iarr,:,:) + log10initval
      xx_genarr2d(iarr,:,:) = EXP(ln10 * xx_genarr2d(iarr,:,:)) 
      fld = xx_genarr2d(iarr,:,:) * xx_genarr2d_mask(iarr,:,:)
    else
      fld = fld + xx_genarr2d(iarr,:,:) * xx_genarr2d_mask(iarr,:,:)
    end if       

  end subroutine ctrl_map_genarr2d

  subroutine ctrl_map_ini_genarr3d

    implicit none

    integer(i4b)        :: ctrl_index
    integer(i4b)        :: igen_temp_c, igen_age_c

    xx_genarr3d_vars            = XX_GENARR3D_VARS_ARR
    xx_genarr3d_preproc         = XX_GENARR3D_PREPROC_ARR
    xx_genarr3d_bounds          = XX_GENARR3D_BOUNDS_ARR
    xx_genarr3d_log10initval    = XX_GENARR3D_LOG10INITVAL_ARR
    xx_genarr3d_weight          = XX_GENARR3D_WEIGHT_ARR

    igen_temp_c = 0
    igen_age_c = 0

    do ctrl_index = 1, NUM_CTRL_GENARR3D
      
      if (trim(adjustl(xx_genarr3d_vars(ctrl_index))) .EQ. 'xx_temp_c') then
        igen_temp_c = ctrl_index
#if (((ANF_DAT == 1) && TEMP_INIT==5) || (ANF_DAT > 2))
        errormsg = ' >>> ctrl_map_ini_genarr3d: ' &
          //'Initial ice temperature as a control param is only compatible with ' &
          //'ANF_DAT==1 and TEMP_INIT<=5 or ANF_DAT==2!'
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

    if (igen_temp_c .GT. 0) then
      call ctrl_map_genarr3d(temp_c, igen_temp_c)
    end if
    if (igen_age_c .GT. 0) then
      call ctrl_map_genarr3d(age_c, igen_age_c)
    end if
    
  end subroutine ctrl_map_ini_genarr3d

  subroutine ctrl_map_genarr3d(fld, iarr)

    implicit none  

    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)          :: fld, wgenarr3d
    integer(i4b)                                        :: iarr, kc, j, i, k3
    integer(i4b)                                        :: iLow, iHigh, jLow, jHigh, kcLow, kcHigh
    logical                                             :: doread, dobounds
    logical                                             :: dosmooth, doscaling, dolog10ctrl
    real(dp)                                            :: log10initval, ln10
    character(128), dimension(NUMCTRLPROCARR3D)            :: preprocs

    ln10 = log(10.0)

    doread      = .FALSE.
    dosmooth    = .FALSE.
    doscaling   = .FALSE.
    dolog10ctrl = .FALSE.
    dobounds    = .FALSE.

    read (unit=xx_genarr3d_preproc(iarr),fmt=*) preprocs

    do k3 = 1, NUMCTRLPROCARR3D

      if (preprocs(k3) .EQ. 'smooth') then
        dosmooth = .TRUE.
        errormsg = ' >>> ctrl_map_genarr3d: ' &
          //'Weaver and Courtier like smoothing is not yet implemented!'
        call error(errormsg)
      end if
      if (preprocs(k3) .EQ. 'scaling') then
        doscaling = .TRUE.
        do i = 0, IMAX
          do j = 0, JMAX
            do kc = 0, KCMAX
              read(unit=xx_genarr3d_weight(iarr), fmt=*) wgenarr3d(kc,j,i)
            end do
          end do
        end do
      else
        wgenarr3d = 1.0
      end if
      if (preprocs(k3) .EQ. 'read') then
        doread = .TRUE.
        errormsg = ' >>> ctrl_map_genarr3d: ' &
          //'Reading xx_gen from a file is yet to be implemented!'
        call error(errormsg)
      end if
      if (preprocs(k3) .EQ. 'log10ctrl') then
        dolog10ctrl = .TRUE.
        log10initval = xx_genarr3d_log10initval(iarr)
      end if
      if (preprocs(k3) .EQ. 'bounds') then
        dobounds = .TRUE.
      end if
    end do

    if (dobounds .AND. trim(adjustl(xx_genarr3d_bounds(iarr))) .NE. ' ') then
      read (unit=xx_genarr3d_bounds(iarr),fmt=*) iLow, iHigh, jLow, jHigh, kcLow, kcHigh
    else
      iLow   = 0
      iHigh  = IMAX
      jLow   = 0
      jHigh  = JMAX
      kcLow  = 0
      kcHigh = KCMAX
    end if

    do i = 0, IMAX
      do j = 0, JMAX
        do kc = 0, KCMAX
        if ((i .GE. iLow) .AND. (i .LE. iHigh) .AND. (j .GE. jLow) .AND. (j .LE. jHigh) &
           .AND. (kc .GE. kcLow) .AND. (kc .LE. kcHigh)) then
          xx_genarr3d_mask(iarr,kc,j,i) = 1.0
        else
          xx_genarr3d_mask(iarr,kc,j,i) = 0.0
        end if
        end do
      end do
    end do

    xx_genarr3d(iarr,:,:,:) = xx_genarr3d(iarr,:,:,:) / sqrt(wgenarr3d(:,:,:))

    if (dolog10ctrl) then  
      xx_genarr3d(iarr,:,:,:) = xx_genarr3d(iarr,:,:,:) + log10initval
      xx_genarr3d(iarr,:,:,:) = EXP(ln10 * xx_genarr3d(iarr,:,:,:)) 
      fld = xx_genarr3d(iarr,:,:,:) * xx_genarr3d_mask(iarr,:,:,:)
    else
      fld = fld + xx_genarr3d(iarr,:,:,:) * xx_genarr3d_mask(iarr,:,:,:)
    end if       

  end subroutine ctrl_map_genarr3d

end module ctrl_map_genarr_m 
#endif /* ALLOW_GENCTRL */
