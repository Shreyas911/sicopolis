!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  o p e n a d _ m
!
!> @file
!!
!! A catch-all module for openad-related subroutines. 
!!
!! @section Copyright
!!
!! Copyright 2017-2021 Liz Curry-Logan, Sri Hari Krishna Narayanan,
!!                     Patrick Heimbach, Ralf Greve
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
!> Module for all openad-related subroutines 
!<------------------------------------------------------------------------------
module openad_m

  implicit none

  private
  public :: adjoint_master

#if   (defined (ALLOW_GRDCHK))
  public :: grdchk_main
#elif (defined (ALLOW_OPENAD))
  private :: direct_substitution
  private :: print_output
#endif

contains

!-------------------------------------------------------------------------------
!> Adjoint master is the main tool by which sicopolis.F90 invokes the adjoint
!! code. Its job is to figure out what mode of the adjoint code is being invoked
!! and run the appropriate subroutine. 
!<------------------------------------------------------------------------------
  subroutine adjoint_master

  implicit none

#ifdef ALLOW_OPENAD
  integer           :: mode
  character(len=32) :: arg
  logical           :: ISPLAIN, ISTAPE, ISADJOINT
#endif

#ifdef ALLOW_GRDCHK
    call grdchk_main
#endif

#ifdef ALLOW_OPENAD
    call get_ad_mode(mode,arg,ISPLAIN, ISTAPE, ISADJOINT)
    call direct_substitution(ISPLAIN, ISTAPE, ISADJOINT)
#endif

  end subroutine adjoint_master

#ifdef ALLOW_GRDCHK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Subroutine :  g r d c h k _ m a i n
!   Purpose    :  Gradient check top level routine
!                 Compares the gradients calculated by the adjoint model
!                 to the gradients calculated by finite difference
!                 approximations
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine grdchk_main
   
   use sico_types_m
   use sico_variables_m
   use sico_vars_m
   
   use sico_init_m
   use sico_main_loop_m
   use sico_end_m
   
   use ctrl_m, only: ctrl_init, cost_independent_init, cost_dependent_init, &
                     cost_final
   
   implicit none
   
   integer(i4b)       :: ndat2d, ndat3d
   integer(i4b)       :: n_output
   real(dp)           :: delta_ts, glac_index
   real(dp)           :: mean_accum
   real(dp)           :: dtime, dtime_temp, dtime_wss, &
                                      dtime_out, dtime_ser
   real(dp)           :: time, time_init, time_end, time_output(100)
   real(dp)           :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
   real(dp)           :: z_sl, dzsl_dtau, z_mar
   character(len=100) :: runname
   
   !-------- Variable declarations needed for this routine specifically
   real(dp)                          :: orig_val, perturb_val = 5e-2
   real(dp),     dimension(3)        :: fc_collected
   real(dp),     dimension(3)        :: direction
   real(dp)                          :: gfd0,gfd, perturbation
   integer(i4b), parameter           :: points = 5
   integer(i4b), dimension(points)   :: ipoints, jpoints
   integer(i4b)                      :: i, j, p, d
   character(len=100)                :: fname
   
   !-------- This array holds the direction of perturbations to follow:
   direction(1) = 0
   direction(2) = 1
   direction(3) = -1

#if (!defined(GRL) && !defined(ANT))
   print *, ">>> Adjoint only available for GRL and ANT right now; kill code." 
#endif

   !-------- Test points along spines of the ice sheets
   do p = 1, points
#if (defined(GRL))
      ipoints(p) = int(real(IMAX/2))
      jpoints(p) = int(real(JMAX/5)) + (p-1) * points
#elif (defined(ANT))
      ipoints(p) = int(real(IMAX/3)) + int(real((.85-.33)*IMAX/points)) * (p - 1) 
      jpoints(p) = int(real(JMAX/2)) 
#endif
   end do

   !-------- Initialize output files 
   open(99, file='GradientVals_'//trim(RUNNAME)//'.dat',&
       form="FORMATTED", status="REPLACE")
   open(98, file='CostVals_'//trim(RUNNAME)//'.dat',&
       form="FORMATTED", status="REPLACE")
   
   !-------- Loop over points
   do p = 1, points
     i = ipoints(p)
     j = jpoints(p)

          !-------- Loop over perturbation direction (0, +, -)
          do d = 1, 3 

          !-------- Let yourself know where you are:
          print *, ' point (p, i, j), direction (d) [ ', p , ', ', i, ', ', j, ', ', d, ' ] '

          !-------- One complete forward run 
            call deldirs
            call cost_dependent_init
        
            call sico_init(delta_ts, glac_index, &
                 mean_accum, &
                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 z_sl, dzsl_dtau, z_mar, &
                 ndat2d, ndat3d, n_output, &
                 runname)

            perturbation = 1 + direction(d) * perturb_val 


          !-------- Controls to be perturbed (add your own here and below in
          !         subroutine print_output()
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir) 
 
            ! -- H_c
            orig_val = H_c(j,i)
            H_c(j,i) = orig_val * perturbation 

            ! -- mean annual temp 
            !orig_val = temp_ma_present(j,i)
            !temp_ma_present(j,i) = orig_val * perturbation

            ! -- mean annual temp 
            !orig_val = q_geo(j,i)
            !q_geo(j,i) = orig_val * perturbation 

            ! -- mean annual temp 
            !orig_val = vx_c(24,j,i)
            !vx_c(24,j,i) = orig_val * perturbation

            ! -- mean annual temp 
            !orig_val = c_slide(j,i)
            !c_slide(j,i) = orig_val * perturbation

            ! -- mean annual temp 
            !orig_val = c_drag(j,i)
            !c_drag(j,i) = orig_val * perturbation

            ! -- precip_present
            !orig_val = precip_present(j,i,1)
            !precip_present(j,i,1) = orig_val * perturbation

            ! -- experimental field acc_fact
            !orig_val = acc_fact(j,i)
            !acc_fact(j,i) = orig_val * perturbation

            ! -- sanity check
            write(6,fmt='(a,f40.20)') "orig_val = ", orig_val
            write(6,fmt='(a,f40.20)') "pert_val = ", orig_val*perturbation

            call sico_main_loop(delta_ts, glac_index, &
                 mean_accum, &
                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 z_sl, dzsl_dtau, z_mar, &
                 ndat2d, ndat3d, n_output, &
                 runname)
          
            call cost_final(runname)
            call sico_end
       
            ! store cost
            fc_collected(d) = fc
        
            ! --------- calculate simple 2-sided finite difference due to
            !           perturbation: fc(+) - fc(-) / 2*espsilon
            gfd = (fc_collected(2) - fc_collected(3))/(2.d0 * perturb_val * orig_val)
          
          end do ! (close perturb loop)

          ! -- sanity check
          write(6, fmt='(a,f40.20)')   "Finite difference is = ", gfd
        
          ! --------- write these values to output file
          write(99, fmt='(f40.20)') gfd
          write(98, fmt='(f40.20)') fc_collected(1)
          write(98, fmt='(f40.20)') fc_collected(2)
          write(98, fmt='(f40.20)') fc_collected(3)
          write(98, fmt='(a)') '----------------------------------'
    
   end do ! (close loop over points)
  
   close(unit=99)
   close(unit=98)
   
   end subroutine grdchk_main
#endif /* Only ALLOW_GRDCHK */

!!-------------------------------------------------------------------------------
!!> Checks to see if output dir exists. If so, deletes it.
!!<------------------------------------------------------------------------------
  subroutine deldirs
  
  implicit none
  
  character(len=256) :: shell_command
  
  !-------- deleting directories
  
    shell_command = 'if [ -d'
    shell_command = trim(shell_command)//' '//OUT_PATH
    shell_command = trim(shell_command)//' '//'] ; then rm -rf'
    shell_command = trim(shell_command)//' '//OUT_PATH
    shell_command = trim(shell_command)//' '//'; fi'
    
    call system(trim(shell_command))
  
  end subroutine deldirs

#ifdef ALLOW_OPENAD
!-------------------------------------------------------------------------------
!> This is the top-level routine for the adjoint mode sweep. It takes arguments
!! from the executable {-p, -t, -a} supplied at run time for the plain, tape, and
!! adjoint modes of execution. 
!<------------------------------------------------------------------------------
  subroutine direct_substitution(ISPLAIN, ISTAPE, ISADJOINT)

  use OAD_rev
  use OAD_tape
  use sico_variables_m
#if (defined(GRL) && DISC>0)
  use discharge_workers_m
#endif
  use ice_material_properties_m
  use enth_temp_omega_m
  use sico_init_m
  use sico_end_m

  logical                                    :: ISPLAIN, ISTAPE, ISADJOINT
  integer(i4b)                               :: ndat2d, ndat3d
  integer(i4b)                               :: n_output
  real(dp)                                   :: delta_ts, glac_index
  real(dp)                                   :: mean_accum
  real(dp)                                   :: dtime, dtime_temp, &
                                                dtime_wss, dtime_out, dtime_ser
  real(dp)                                   :: time, time_init, time_end
  real(dp), dimension(100)                   :: time_output
  real(dp)                                   :: dxi, deta, dzeta_c, &
                                                dzeta_t, dzeta_r
  real(dp)                                   :: z_sl, dzsl_dtau, z_mar
  character(len=100)                         :: runname

  our_rev_mode%arg_store  = .false.
  our_rev_mode%arg_restore= .false.

  our_rev_mode%res_store  = .false.
  our_rev_mode%res_restore= .false.

  our_rev_mode%plain      = .false.
  our_rev_mode%tape       = .false.
  our_rev_mode%adjoint    = .false.

 call oad_independent_init

 call sico_init(delta_ts, glac_index, &
      mean_accum, &
      dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
      time, time_init, time_end, time_output, &
      dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
      z_sl, dzsl_dtau, z_mar, &
      ndat2d, ndat3d, n_output, &
      runname)

 call var_transfer&
      (maske, maske_old, maske_neu, &
      n_cts, n_cts_neu, kc_cts, kc_cts_neu, &
      flag_calc_temp, &
      flag_inner_point, &
      flag_grounding_line_1,flag_grounding_line_2,flag_calving_front_1,&
      flag_calving_front_2,flag_shelfy_stream_x,flag_shelfy_stream_y,&
      flag_shelfy_stream,xi,eta,zeta_c,zeta_t,zeta_r,aa,flag_aa_nonzero,&
      ea,eaz_c,eaz_c_quotient,lambda,phi,area,sq_g11_g,sq_g22_g,&
      insq_g11_g,insq_g22_g,sq_g11_sgx,sq_g11_sgy,sq_g22_sgx,sq_g22_sgy,&
      insq_g11_sgx,insq_g22_sgy,zs,zm,zb,zl,zl0,wss,flex_rig_lith,&
      time_lag_asth,H_c,H_t,dzs_dxi,dzm_dxi,dzb_dxi,dH_c_dxi,dH_t_dxi,&
      dzs_deta,dzm_deta,dzb_deta,dH_c_deta,dH_t_deta,dzs_dxi_g,&
      dzm_dxi_g,dzb_dxi_g,dH_c_dxi_g,dH_t_dxi_g,dzs_deta_g,dzm_deta_g,&
      dzb_deta_g,dH_c_deta_g,dH_t_deta_g,dzs_dtau,dzm_dtau,dzb_dtau,&
      dzl_dtau,dH_c_dtau,dH_t_dtau,p_weert,q_weert,p_weert_inv,&
      c_slide,d_help_b,c_drag,p_b_w,&
      vx_b,vy_b,&
      vx_m,vy_m,vx_m_sia,vy_m_sia,vx_m_ssa,vy_m_ssa,&
      ratio_sl_x,ratio_sl_y,ratio_sl,&
      vx_b_g,vy_b_g,vz_b,vz_m,vx_s_g,vy_s_g,vz_s,&
      flui_ave_sia,h_diff,qx,qy,q_gl_g,q_geo,temp_b,temph_b,Q_bm,Q_b_apl,&
      Q_tld,Q_b_tot,H_w,&
      accum,runoff,runoff_apl,as_perp,temp_maat,temp_s,am_perp,&
      am_perp_st,zs_neu,zm_neu,zb_neu,zl_neu,H_c_neu,H_t_neu,zs_ref,&
      accum_present,precip_ma_present,precip_ma_lgm_anom,&
      temp_ma_present,temp_mj_present,temp_ma_lgm_anom,temp_mj_lgm_anom,&
      dist_dxdy,acc_fact,precip_present,precip_lgm_anom,gamma_precip_lgm_anom,&
      temp_mm_present,temp_mm_lgm_anom,d_help_c,vx_c,vy_c,vz_c,temp_c,&
      temp_c_neu,temp_c_m,age_c,age_c_neu,txz_c,tyz_c,sigma_c,enh_c,&
      de_ssa,vis_int_g,vx_g,vy_g,d_help_t,vx_t,vy_t,vz_t,omega_t,&
      omega_t_neu,temp_t_m,age_t,age_t_neu,txz_t,tyz_t,sigma_t,enh_t,&
      temp_r,temp_r_neu,enth_c,enth_c_neu,omega_c,omega_c_neu,enth_t,&
      enth_t_neu,dxx_c,dyy_c,dxy_c,dxz_c,dyz_c,de_c,lambda_shear_c,&
      dxx_t,dyy_t,dxy_t,dxz_t,dyz_t,de_t,lambda_shear_t,RHO,RHO_W,&
      RHO_SW,L,G,NUE,BETA,DELTA_TM_SW,OMEGA_MAX,H_R,RHO_C_R,KAPPA,KAPPA_R,&
      RHO_A,R_T,R_T_IMQ,&
      R,A,B,F_INV,LAMBDA0,PHI0,&
      S_STAT_0,BETA1_0,BETA1_LT_0,&
      BETA1_HT_0,BETA2_0,BETA2_LT_0,BETA2_HT_0,PHI_SEP_0,PMAX_0,MU_0,&
      RF,RF_IMQ,KAPPA_IMQ,C_IMQ,year_zero,&
      forcing_flag,n_core,lambda_core,phi_core,x_core,&
      y_core,grip_time_min,grip_time_stp,grip_time_max,ndata_grip,&
      griptemp,gi_time_min,gi_time_stp,gi_time_max,ndata_gi,&
      glacial_index,specmap_time_min,specmap_time_stp,specmap_time_max,&
      ndata_specmap,specmap_zsl,time_target_topo_init,&
      time_target_topo_final,maske_target,zs_target,zb_target,&
      zl_target,H_target,maske_maxextent,mask_ablation_type,ncid_temp_precip,&
      ndata_temp_precip,temp_precip_time_min,temp_precip_time_stp,&
      temp_precip_time_max,temp_precip_time,kei,kei_r_max,&
      kei_r_incr,fc,objf_test,&
      mult_test,target_topo_tau_0,&
      c_int_table,c_int_inv_table,C,n_temp_min,n_temp_max,n_enth_min,&
      n_enth_max,L_inv,L_eto,&
      RHO_I_IMQ,RHO_C_IMQ,KAPPA_C_IMQ,C_C_IMQ,&
      as_perp_apl,mb_source_apl,accum_apl,&
      calving,calving_apl,ncid_ser,ncid_core,n_data_kei,&
#if (defined(INITMIP_SMB_ANOM_FILE))
      as_anom_initmip,ab_anom_initmip,&
#endif
      n_slide_region,&
#if (DISC==2)
      glann_time_min,glann_time_stp,glann_time_max,&
      ndata_glann,dT_glann_CLIMBER,&
#endif
#if (defined(GRL) && DISC>0)
      disc_DW,n_discharge_call_DW,iter_mar_coa_DW,c_dis_0_DW,s_dis_DW,&
      c_dis_fac_DW,T_sub_PD_DW,alpha_sub_DW,alpha_o_DW,m_H_DW,&
      m_D_DW,r_mar_eff_DW,T_sea_freeze_DW,c_dis_DW,&
#elif (defined(ANT))
      n_bm_region,&
#endif
      et,melt,melt_star,rainfall,snowfall,&
#if (defined(AGE_COST))
      age_data,age_unc, &
#endif 
      n_temp_min_IMQ,n_temp_max_IMQ,&
      smb_corr,smb_corr_prescribed,&
      flag_grounded_front_a_1,flag_grounded_front_a_2,&
      flag_grounded_front_b_1,flag_grounded_front_b_2,&
      ij2n,n2i,n2j,&
      Q_w,Q_w_x,Q_w_y,&
      vis_ave_g)

      our_rev_mode%plain=ISPLAIN
      our_rev_mode%tape=ISTAPE

      if(ISTAPE.eqv..true.) call oad_tape_init()
      call sicopolis_openad(delta_ts, glac_index, &
      mean_accum, &
      dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
      time, time_init, time_end, time_output, &
      dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
      z_sl, dzsl_dtau, z_mar, &
      ndat2d, ndat3d, n_output, &
      runname)

      if(ISADJOINT) then
        our_rev_mode%tape=.FALSE.
        our_rev_mode%adjoint=.TRUE.
        call sicopolis_openad(delta_ts, glac_index, &
        mean_accum, &
        dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
        time, time_init, time_end, time_output, &
        dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
        z_sl, dzsl_dtau, z_mar, &
        ndat2d, ndat3d, n_output, &
        runname)
        call print_output(runname)
        call oad_tape_delete()
      else 
        if(ISTAPE.eqv..true.) call oad_tape_delete()
      end if

      call sico_end

  end subroutine direct_substitution

!-------------------------------------------------------------------------------
!> This is where the final cost fc is initialized. See ctrl_m.F90 (subroutine
!! final_cost) to see its structure. 
!<------------------------------------------------------------------------------
  subroutine oad_independent_init()

  use OAD_active
  use oad_sico_variables_m

  implicit none

  fc%d = 1.0_dp

  end subroutine oad_independent_init

!-------------------------------------------------------------------------------
!> Take arguments passed to executable called driver in OAD mode.
!! Only options are (should be)
!! "-h help"
!! "-p PLAIN MODE"
!! "-t TAPE MODE"
!! "-a ADJOINT MODE" 
!!
!<------------------------------------------------------------------------------
  subroutine get_ad_mode(mode,arg,ISPLAIN, ISTAPE, ISADJOINT)

  implicit none

  integer           :: mode
  character(len=32) :: arg

  logical :: ISPLAIN, ISTAPE, ISADJOINT

  ISPLAIN   = .false.
  ISTAPE    = .false.
  ISADJOINT = .false.

  do mode=1,command_argument_count()
    ! grabbing the input to driver
    call get_command_argument(mode,arg)

    ! going through possible arguments to driver
    select case(adjustl(arg))
      case("-p")
         ISPLAIN = .TRUE.
      case("-t")
         ISTAPE = .TRUE.
      case("-a")
         ISADJOINT = .TRUE.
      case("-h")
         print *, "Command line options:"
         print *, "-h help"
         print *, "-p PLAIN MODE"
         print *, "-t TAPE MODE"
         print *, "-a ADJOINT MODE"
         stop
      case default
         print *, "Unknown command line option ", arg
         print *, "Command line options are:"
         print *, "-h help"
         print *, "-p PLAIN MODE"
         print *, "-t TAPE MODE"
         print *, "-a ADJOINT MODE"
         stop
    end select
  end do

  if (    (ISPLAIN   .NEQV. .TRUE.) &
    .AND. (ISTAPE    .NEQV. .TRUE.) &
    .AND. (ISADJOINT .NEQV. .TRUE.)) then
      print *, " No valid option specified."
      print *, "Command line options are:"
      print *, "-h help"
      print *, "-p PLAIN MODE"
      print *, "-t TAPE MODE"
      print *, "-a ADJOINT MODE"
      stop
  end if

  if(ISADJOINT .EQV. .TRUE.) then
    ISTAPE = .TRUE.
  end if

  if((ISPLAIN .EQV. .TRUE.) .AND. (ISTAPE .EQV. .TRUE.)) then
    print *, " Cannot specify both -p along with any other option"
    stop
  end if

  end subroutine get_ad_mode 

!-------------------------------------------------------------------------------
!> This is an absolutely essential routine that takes all of the variables
!! foobar used in the forward (tape) sweep and assigns them to a_foobar
!! variables *needed* to perform the reverse (adjoint) sweep. If your
!! compilation FAILS with the error message complaining that something was not
!! ACTIVE or was ACTIVE and needed not be, the solution is to be found inside
!! this routine, by taking variable foobar and making it either active by
!! assigning:
!! foobar%v = a_foobar (to make foobar ACTIVE) or
!! foobar   = a_foobar (to make foobar NOT ACTIVE). But that's not all!
!! If your regression script suddenly is doing quite poorly (FDs vs. ADs are not
!! within tolerance anymore) is it likely that some new variable as been
!! introduced in the trunk and needs to be *transfered* here. In other words, in
!! the commit history of sicopolis, someone introduced new_foobar, and
!! new_foobar is completely missing in the below *and* in the calling (in
!! subroutine direct_substitution()) of var_transfer above. This needs to be
!! amended before proceeding. If you feel unsure, email
!! liz.curry.logan@gmail.com. 
!<------------------------------------------------------------------------------
  subroutine var_transfer&
         (a_maske, a_maske_old, a_maske_neu, &
         a_n_cts, a_n_cts_neu, a_kc_cts, a_kc_cts_neu, &
         a_flag_calc_temp, &
         a_flag_inner_point, &
         a_flag_grounding_line_1,a_flag_grounding_line_2,a_flag_calving_front_1,&
         a_flag_calving_front_2,a_flag_shelfy_stream_x,a_flag_shelfy_stream_y,&
         a_flag_shelfy_stream,a_xi,a_eta,a_zeta_c,a_zeta_t,a_zeta_r,a_aa,a_flag_aa_nonzero,&
         a_ea,a_eaz_c,a_eaz_c_quotient,a_lambda,a_phi,a_area,a_sq_g11_g,a_sq_g22_g,&
         a_insq_g11_g,a_insq_g22_g,a_sq_g11_sgx,a_sq_g11_sgy,a_sq_g22_sgx,a_sq_g22_sgy,&
         a_insq_g11_sgx,a_insq_g22_sgy,a_zs,a_zm,a_zb,a_zl,a_zl0,a_wss,a_flex_rig_lith,&
         a_time_lag_asth,a_H_c,a_H_t,a_dzs_dxi,a_dzm_dxi,a_dzb_dxi,a_dH_c_dxi,a_dH_t_dxi,&
         a_dzs_deta,a_dzm_deta,a_dzb_deta,a_dH_c_deta,a_dH_t_deta,a_dzs_dxi_g,&
         a_dzm_dxi_g,a_dzb_dxi_g,a_dH_c_dxi_g,a_dH_t_dxi_g,a_dzs_deta_g,a_dzm_deta_g,&
         a_dzb_deta_g,a_dH_c_deta_g,a_dH_t_deta_g,a_dzs_dtau,a_dzm_dtau,a_dzb_dtau,&
         a_dzl_dtau,a_dH_c_dtau,a_dH_t_dtau,a_p_weert,a_q_weert,a_p_weert_inv,&
         a_c_slide,a_d_help_b,a_c_drag,a_p_b_w,&
         a_vx_b,a_vy_b,&
         a_vx_m,a_vy_m,a_vx_m_sia,a_vy_m_sia,a_vx_m_ssa,a_vy_m_ssa,&
         a_ratio_sl_x,a_ratio_sl_y,a_ratio_sl,&
         a_vx_b_g,a_vy_b_g,a_vz_b,a_vz_m,a_vx_s_g,a_vy_s_g,a_vz_s,&
         a_flui_ave_sia,a_h_diff,a_qx,a_qy,a_q_gl_g,a_q_geo,a_temp_b,a_temph_b,a_Q_bm,a_Q_b_apl,&
         a_Q_tld,a_Q_b_tot,a_H_w,&
         a_accum,a_runoff,a_runoff_apl,a_as_perp,a_temp_maat,a_temp_s,a_am_perp,&
         a_am_perp_st,a_zs_neu,a_zm_neu,a_zb_neu,a_zl_neu,a_H_c_neu,a_H_t_neu,a_zs_ref,&
         a_accum_present,a_precip_ma_present,a_precip_ma_lgm_anom,&
         a_temp_ma_present,a_temp_mj_present,a_temp_ma_lgm_anom,a_temp_mj_lgm_anom,&
         a_dist_dxdy,a_acc_fact,a_precip_present,a_precip_lgm_anom,a_gamma_precip_lgm_anom,&
         a_temp_mm_present,a_temp_mm_lgm_anom,a_d_help_c,a_vx_c,a_vy_c,a_vz_c,a_temp_c,&
         a_temp_c_neu,a_temp_c_m,a_age_c,a_age_c_neu,a_txz_c,a_tyz_c,a_sigma_c,a_enh_c,&
         a_de_ssa,a_vis_int_g,a_vx_g,a_vy_g,a_d_help_t,a_vx_t,a_vy_t,a_vz_t,a_omega_t,&
         a_omega_t_neu,a_temp_t_m,a_age_t,a_age_t_neu,a_txz_t,a_tyz_t,a_sigma_t,a_enh_t,&
         a_temp_r,a_temp_r_neu,a_enth_c,a_enth_c_neu,a_omega_c,a_omega_c_neu,a_enth_t,&
         a_enth_t_neu,a_dxx_c,a_dyy_c,a_dxy_c,a_dxz_c,a_dyz_c,a_de_c,a_lambda_shear_c,&
         a_dxx_t,a_dyy_t,a_dxy_t,a_dxz_t,a_dyz_t,a_de_t,a_lambda_shear_t,a_RHO,a_RHO_W,&
         a_RHO_SW,a_L,a_G,a_NUE,a_BETA,a_DELTA_TM_SW,a_OMEGA_MAX,a_H_R,a_RHO_C_R,a_KAPPA,a_KAPPA_R,&
         a_RHO_A,a_R_T,a_R_T_IMQ,&
         a_R,a_A,a_B,a_F_INV,a_LAMBDA0,a_PHI0,&
         a_S_STAT_0,a_BETA1_0,a_BETA1_LT_0,&
         a_BETA1_HT_0,a_BETA2_0,a_BETA2_LT_0,a_BETA2_HT_0,a_PHI_SEP_0,a_PMAX_0,a_MU_0,&
         a_RF,a_RF_IMQ,a_KAPPA_IMQ,a_C_IMQ,a_year_zero,&
         a_forcing_flag,a_n_core,a_lambda_core,a_phi_core,a_x_core,&
         a_y_core,a_grip_time_min,a_grip_time_stp,a_grip_time_max,a_ndata_grip,&
         a_griptemp,a_gi_time_min,a_gi_time_stp,a_gi_time_max,a_ndata_gi,&
         a_glacial_index,a_specmap_time_min,a_specmap_time_stp,a_specmap_time_max,&
         a_ndata_specmap,a_specmap_zsl,a_time_target_topo_init,&
         a_time_target_topo_final,a_maske_target,a_zs_target,a_zb_target,&
         a_zl_target,a_H_target,a_maske_maxextent,a_mask_ablation_type,a_ncid_temp_precip,&
         a_ndata_temp_precip,a_temp_precip_time_min,a_temp_precip_time_stp,&
         a_temp_precip_time_max,a_temp_precip_time,a_kei,a_kei_r_max,&
         a_kei_r_incr,a_fc,a_objf_test,&
         a_mult_test,a_target_topo_tau_0,&
         a_c_int_table,a_c_int_inv_table,a_c,a_n_temp_min,a_n_temp_max,a_n_enth_min,&
         a_n_enth_max,a_L_inv,a_L_eto,&
         a_RHO_I_IMQ,a_RHO_C_IMQ,a_KAPPA_C_IMQ,a_C_C_IMQ, &
         a_as_perp_apl,a_mb_source_apl,a_accum_apl, &
         a_calving,a_calving_apl,a_ncid_ser,a_ncid_core,a_n_data_kei,&
#if (defined(INITMIP_SMB_ANOM_FILE))
         a_as_anom_initmip,a_ab_anom_initmip,&
#endif
         a_n_slide_region,&
#if (DISC==2)
         a_glann_time_min,a_glann_time_stp,a_glann_time_max,&
         a_ndata_glann,a_dT_glann_CLIMBER,&
#endif
#if (defined(GRL) && DISC>0)
         a_disc_DW,a_n_discharge_call_DW,a_iter_mar_coa_DW,a_c_dis_0_DW,a_s_dis_DW,&
         a_c_dis_fac_DW,a_T_sub_PD_DW,a_alpha_sub_DW,a_alpha_o_DW,a_m_H_DW,&
         a_m_D_DW,a_r_mar_eff_DW,a_T_sea_freeze_DW,a_c_dis_DW,&
#elif (defined(ANT))
         a_n_bm_region,&
#endif
         a_et,a_melt,a_melt_star,a_rainfall,a_snowfall,&
#if (defined(AGE_COST))
         a_age_data,a_age_unc, &
#endif 
         a_n_temp_min_IMQ,a_n_temp_max_IMQ,&
         a_smb_corr,a_smb_corr_prescribed,&
         a_flag_grounded_front_a_1,a_flag_grounded_front_a_2,&
         a_flag_grounded_front_b_1,a_flag_grounded_front_b_2,&
         a_ij2n,a_n2i,a_n2j,&
         a_Q_w,a_Q_w_x,a_Q_w_y,&
         a_vis_ave_g)
  
    use OAD_active
    use oad_sico_variables_m
    use oad_sico_vars_m
    use oad_enth_temp_omega_m
    use oad_ice_material_properties_m
#if (defined(GRL) && DISC>0)
    use oad_discharge_workers_m
#endif
  
    implicit none
  
    real(dp)                                           :: a_A
    real(dp)                                           :: a_aa
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_acc_fact
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_accum
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_accum_present
#if (defined(INITMIP_SMB_ANOM_FILE))
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_ab_anom_initmip
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_as_anom_initmip
#endif
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_age_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_age_c_neu
#if (defined(AGE_COST))
#if (CALCMOD==1)
    real(dp), dimension(0:KTMAX+KCMAX+1,0:JMAX,0:IMAX) :: a_age_data
    real(dp), dimension(0:KTMAX+KCMAX+1,0:JMAX,0:IMAX) :: a_age_unc
#else
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_age_data
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_age_unc
#endif
#endif /* No age cost used */
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_age_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_age_t_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_area
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_am_perp
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_am_perp_st
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_as_perp
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_as_perp_apl
    real(dp)                                           :: a_B
    real(dp)                                           :: a_BETA
    real(dp)                                           :: a_BETA1_0
    real(dp)                                           :: a_BETA1_HT_0
    real(dp)                                           :: a_BETA1_LT_0
    real(dp)                                           :: a_BETA2_0
    real(dp)                                           :: a_BETA2_HT_0
    real(dp)                                           :: a_BETA2_LT_0
    real(dp), dimension(-190:10)                       :: a_c
    real(dp)                                           :: a_C_C_IMQ
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_c_drag
    real(dp), dimension(-256:255)                      :: a_C_IMQ
    real(dp), dimension(-256:255)                      :: a_c_int_table
    real(dp), dimension(-524288:524287)                :: a_c_int_inv_table
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_c_slide
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_calving
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_calving_apl
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_d_help_b
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_d_help_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_d_help_t
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_de_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_de_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_de_ssa
    real(dp)                                           :: a_DELTA_TM_SW
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_c_deta
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_c_deta_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_c_dtau
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_c_dxi
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_c_dxi_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_t_deta
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_t_deta_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_t_dtau
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_t_dxi
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dH_t_dxi_g
    real(dp), dimension(-JMAX:JMAX,-IMAX:IMAX)         :: a_dist_dxdy
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_dxx_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_dxy_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_dxz_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_dyy_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_dyz_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_dxx_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_dxy_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_dxz_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_dyy_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_dyz_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzb_deta
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzb_deta_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzb_dtau
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzb_dxi
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzb_dxi_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzl_dtau
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzm_deta
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzm_deta_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzm_dtau
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzm_dxi
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzm_dxi_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzs_deta
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzs_deta_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzs_dtau
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzs_dxi
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_dzs_dxi_g
    real(dp)                                           :: a_ea
    real(dp), dimension(0:KCMAX)                       :: a_eaz_c
    real(dp), dimension(0:KCMAX)                       :: a_eaz_c_quotient
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_enh_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_enh_t
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_enth_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_enth_c_neu
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_enth_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_enth_t_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_et
    real(dp), dimension(0:JMAX)                        :: a_eta
    real(dp)                                           :: a_F_INV
    real(dp)                                           :: a_fc
    logical                                            :: a_flag_aa_nonzero
    logical                                            :: a_flag_calc_temp
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_inner_point
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_grounding_line_1
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_grounding_line_2
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_grounded_front_a_1
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_grounded_front_a_2
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_grounded_front_b_1
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_grounded_front_b_2
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_calving_front_1
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_calving_front_2
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_shelfy_stream_x
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_shelfy_stream_y
    logical, dimension(0:JMAX,0:IMAX)                  :: a_flag_shelfy_stream
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_flex_rig_lith
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_flui_ave_sia
    integer(i4b)                                       :: a_forcing_flag
    real(dp)                                           :: a_G
    real(dp), dimension(0:JMAX,0:IMAX,12)              :: a_gamma_precip_lgm_anom
    integer(i4b)                                       :: a_gi_time_max
    integer(i4b)                                       :: a_gi_time_min
    integer(i4b)                                       :: a_gi_time_stp
#if (DISC==2)
    integer(i4b)                                       :: a_glann_time_min
    integer(i4b)                                       :: a_glann_time_stp
    integer(i4b)                                       :: a_glann_time_max
    integer(i4b)                                       :: a_ndata_glann
    real(dp), dimension(:), allocatable                :: a_dT_glann_CLIMBER
#endif
    integer(i4b)                                       :: a_ndata_gi
    real(dp), dimension(0:a_ndata_gi)                  :: a_glacial_index
    integer(i4b)                                       :: a_grip_time_max
    integer(i4b)                                       :: a_grip_time_min
    integer(i4b)                                       :: a_grip_time_stp
    real(dp), dimension(0:a_ndata_grip)                :: a_griptemp
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_h_diff
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_H_c
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_H_c_neu
    real(dp)                                           :: a_H_R
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_H_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_H_target
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_H_t_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_H_w
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_ij2n
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_insq_g11_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_insq_g22_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_insq_g11_sgx
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_insq_g22_sgy
    real(dp), dimension(-190:10)                       :: a_KAPPA
    real(dp)                                           :: a_KAPPA_C_IMQ
    real(dp), dimension(-256:255)                      :: a_KAPPA_IMQ
    real(dp)                                           :: a_KAPPA_R
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_kc_cts
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_kc_cts_neu
    real(dp), dimension(-10000:10000)                  :: a_kei
    real(dp)                                           :: a_kei_r_max
    real(dp)                                           :: a_kei_r_incr
    real(dp)                                           :: a_L
    real(dp)                                           :: a_L_eto
    real(dp)                                           :: a_L_inv
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_lambda
    real(dp)                                           :: a_LAMBDA0
    real(dp), dimension(a_n_core)                      :: a_lambda_core
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_lambda_shear_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_lambda_shear_t
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_maske
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_mask_ablation_type
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_maske_maxextent
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_maske_neu
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_maske_old
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_n_slide_region
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_maske_target
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_melt
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_melt_star
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_mb_source_apl
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_accum_apl
    real(dp)                                           :: a_MU_0
    real(dp)                                           :: a_mult_test
    integer(i4b), dimension((IMAX+1)*(JMAX+1))         :: a_n2i
    integer(i4b), dimension((IMAX+1)*(JMAX+1))         :: a_n2j
    integer(i4b)                                       :: a_n_core
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_n_cts
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_n_cts_neu
    integer(i4b)                                       :: a_n_enth_min
    integer(i4b)                                       :: a_n_enth_max
    integer(i4b)                                       :: a_n_data_kei
#if (defined(ANT))
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_n_bm_region
#endif
    integer(i4b)                                       :: a_n_temp_min
    integer(i4b)                                       :: a_n_temp_max
    integer(i4b)                                       :: a_n_temp_min_IMQ
    integer(i4b)                                       :: a_n_temp_max_IMQ
    integer(i4b)                                       :: a_ncid_core
    integer(i4b)                                       :: a_ncid_ser
    integer(i4b)                                       :: a_ncid_temp_precip
    integer(i4b)                                       :: a_ndata_grip
    integer(i4b)                                       :: a_ndata_specmap
    integer(i4b)                                       :: a_ndata_temp_precip
    real(dp)                                           :: a_NUE
    real(dp)                                           :: a_objf_test
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_omega_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_omega_c_neu
    real(dp)                                           :: a_OMEGA_MAX
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_omega_t
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_omega_t_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_p_b_w
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_p_weert
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_p_weert_inv
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_phi
    real(dp)                                           :: a_PHI0
    real(dp), dimension(a_n_core)                      :: a_phi_core
    real(dp)                                           :: a_PHI_SEP_0
    real(dp)                                           :: a_PMAX_0
    real(dp), dimension(0:JMAX,0:IMAX,12)              :: a_precip_lgm_anom
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_precip_ma_lgm_anom
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_precip_ma_present
    real(dp), dimension(0:JMAX,0:IMAX,12)              :: a_precip_present
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_b_apl
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_b_tot
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_bm
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_q_imp
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_q_geo
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_q_gl_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_tld
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_w
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_w_x
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_Q_w_y
    integer(i4b), dimension(0:JMAX,0:IMAX)             :: a_q_weert
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_qx
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_qy
    real(dp)                                           :: a_R
    real(dp)                                           :: a_R_T
    real(dp)                                           :: a_R_T_IMQ
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_rainfall
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_ratio_sl_x
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_ratio_sl_y
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_ratio_sl
    real(dp), dimension(-190:10)                       :: a_RF
    real(dp), dimension(-256:255)                      :: a_RF_IMQ
    real(dp)                                           :: a_RHO
    real(dp)                                           :: a_RHO_A
    real(dp)                                           :: a_RHO_C_IMQ
    real(dp)                                           :: a_RHO_C_R
    real(dp)                                           :: a_RHO_I_IMQ
    real(dp)                                           :: a_RHO_SW
    real(dp)                                           :: a_RHO_W
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_runoff
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_runoff_apl
    real(dp)                                           :: a_S_STAT_0
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_sigma_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_sigma_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_smb_corr
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_smb_corr_prescribed
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_snowfall
    integer(i4b)                                       :: a_specmap_time_max
    integer(i4b)                                       :: a_specmap_time_min
    integer(i4b)                                       :: a_specmap_time_stp
    real(dp), dimension(0:a_ndata_specmap)             :: a_specmap_zsl
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_sq_g11_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_sq_g22_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_sq_g11_sgx
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_sq_g11_sgy
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_sq_g22_sgx
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_sq_g22_sgy
    real(dp)                                           :: a_target_topo_tau_0
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_temp_c
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_temp_c_m
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_temp_c_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_b
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_maat
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_s
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temph_b
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_ma_lgm_anom
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_ma_present
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_mj_lgm_anom
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_temp_mj_present
    real(dp), dimension(0:JMAX,0:IMAX,12)              :: a_temp_mm_lgm_anom
    real(dp), dimension(0:JMAX,0:IMAX,12)              :: a_temp_mm_present
    real(dp), dimension(0:a_ndata_temp_precip)         :: a_temp_precip_time
    real(dp)                                           :: a_temp_precip_time_max
    real(dp)                                           :: a_temp_precip_time_min
    real(dp)                                           :: a_temp_precip_time_stp
    real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX)         :: a_temp_r
    real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX)         :: a_temp_r_neu
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_temp_t_m
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_time_lag_asth
    real(dp)                                           :: a_time_target_topo_final
    real(dp)                                           :: a_time_target_topo_init
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_txz_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_txz_t
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_tyz_c
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_tyz_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vis_ave_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vis_int_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_b
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_b_g
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_vx_c
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_m
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_m_sia
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_m_ssa
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vx_s_g
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_vx_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_b
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_b_g
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_vy_c
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_g
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_m
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_m_sia
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_m_ssa
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vy_s_g
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_vy_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vz_b
    real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX)         :: a_vz_c
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vz_m
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_vz_s
    real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX)         :: a_vz_t
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_wss
    real(dp), dimension(a_n_core)                      :: a_x_core
    real(dp), dimension(0:IMAX)                        :: a_xi
    real(dp), dimension(a_n_core)                      :: a_y_core
    real(dp)                                           :: a_year_zero
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zb
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zb_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zb_target
    real(dp), dimension(0:KCMAX)                       :: a_zeta_c
    real(dp), dimension(0:KTMAX)                       :: a_zeta_t
    real(dp), dimension(0:KRMAX)                       :: a_zeta_r
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zl
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zl_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zl_target
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zl0
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zm
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zm_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zs
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zs_neu
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zs_ref
    real(dp), dimension(0:JMAX,0:IMAX)                 :: a_zs_target

#if (defined(GRL) && DISC>0)
    integer(i4b) :: a_disc_DW
    integer(i4b) :: a_n_discharge_call_DW
    integer(i4b) :: a_iter_mar_coa_DW
    real(dp)     :: a_c_dis_0_DW
    real(dp)     :: a_s_dis_DW
    real(dp)     :: a_c_dis_fac_DW
    real(dp)     :: a_T_sub_PD_DW
    real(dp)     :: a_alpha_sub_DW
    real(dp)     :: a_alpha_o_DW
    real(dp)     :: a_m_H_DW
    real(dp)     :: a_m_D_DW
    real(dp)     :: a_r_mar_eff_DW
    real(dp)     :: a_T_sea_freeze_DW
    real(dp), dimension(0:JMAX,0:IMAX) :: a_c_dis_DW
#endif

    A = a_A
    aa = a_aa
#if (defined(INITMIP_BMB_ANOM_FILE))
    ab_anom_initmip = a_ab_anom_initmip
#endif
    accum%v = a_accum
    accum_present = a_accum_present
    acc_fact = a_acc_fact
    age_c%v = a_age_c
    age_c_neu%v = a_age_c_neu
#if (defined(AGE_COST))
    age_data = a_age_data
    age_unc = a_age_unc
#endif 
#if (defined(AGE_COST) && CALCMOD==1) 
    age_t%v = a_age_t
    age_t_neu%v = a_age_t_neu
#else
    age_t = a_age_t
    age_t_neu = a_age_t_neu
#endif
#if (defined(GRL) && DISC>0)
    alpha_o_DW = a_alpha_o_DW     !ns
    alpha_sub_DW = a_alpha_sub_DW !ns
#endif
    am_perp = a_am_perp
    am_perp_st = a_am_perp_st
    area = a_area
    as_perp%v = a_as_perp
    as_perp_apl = a_as_perp_apl
#if (defined(INITMIP_SMB_ANOM_FILE))
    as_anom_initmip = a_as_anom_initmip
#endif
    B = a_B
    BETA = a_BETA
    BETA1_0 = a_BETA1_0
    BETA1_HT_0 = a_BETA1_HT_0
    BETA1_LT_0 = a_BETA1_LT_0
    BETA2_0 = a_BETA2_0
    BETA2_HT_0 = a_BETA2_HT_0
    BETA2_LT_0 = a_BETA2_LT_0
    c = a_c
    calving%v = a_calving
    calving_apl = a_calving_apl
    C_IMQ = a_C_IMQ !ns
    C_C_IMQ = a_C_C_IMQ !ns
    c_int_table = a_c_int_table !ns
    c_int_inv_table = a_c_int_inv_table !ns
#if (defined(GRL) && DISC>0)
    c_dis_DW = a_c_dis_DW !ns
    c_dis_0_DW = a_c_dis_0_DW !ns
#endif
#if (DYNAMICS==2 || MARGIN==3)
    c_drag%v = a_c_drag
#else
    c_drag = a_c_drag
#endif
    c_slide%v = a_c_slide
    DELTA_TM_SW = a_DELTA_TM_SW
#if (ENHMOD == 5 || DYNAMICS==2)
    de_c%v = a_de_c
#else
    de_c = a_de_c
#endif
#if (MARGIN==3 || DYNAMICS==2)
    de_ssa%v = a_de_ssa
#else
    de_ssa = a_de_ssa
#endif
    de_t = a_de_t
    dH_c_deta = a_dH_c_deta
    dH_c_deta_g%v = a_dH_c_deta_g
    dH_c_dtau%v = a_dH_c_dtau
    dH_c_dxi = a_dH_c_dxi
    dH_c_dxi_g%v = a_dH_c_dxi_g
    dH_t_deta = a_dH_t_deta
    dH_t_deta_g%v = a_dH_t_deta_g
    dH_t_dtau%v = a_dH_t_dtau
    dH_t_dxi = a_dH_t_dxi
    dH_t_dxi_g%v = a_dH_t_dxi_g
#if (defined(GRL) && DISC>0)
    disc_DW = a_disc_DW !ns
#endif
    dist_dxdy = a_dist_dxdy
#if (CALCMOD == 5 || ENHMOD == 5 || DYNAMICS==2)
    dxx_c%v = a_dxx_c
    dxy_c%v = a_dxy_c
    dxz_c%v = a_dxz_c
    dyy_c%v = a_dyy_c
    dyz_c%v = a_dyz_c
#else
    dxx_c = a_dxx_c
    dxy_c = a_dxy_c
    dxz_c = a_dxz_c
    dyy_c = a_dyy_c
    dyz_c = a_dyz_c
#endif
    dxx_t = a_dxx_t
    dxy_t = a_dxy_t
    dxz_t = a_dxz_t
    dyy_t = a_dyy_t
    dyz_t = a_dyz_t
    dzb_deta = a_dzb_deta
    dzb_deta_g%v = a_dzb_deta_g
    dzb_dtau%v = a_dzb_dtau
    dzb_dxi = a_dzb_dxi
    dzb_dxi_g%v = a_dzb_dxi_g
#if (REBOUND >= 1)
    dzl_dtau%v = a_dzl_dtau
#else
    dzl_dtau = a_dzl_dtau
#endif
    dzm_deta = a_dzm_deta
    dzm_deta_g%v = a_dzm_deta_g
    dzm_dtau%v = a_dzm_dtau
    dzm_dxi = a_dzm_dxi
    dzm_dxi_g%v = a_dzm_dxi_g
    dzs_deta%v = a_dzs_deta
    dzs_deta_g%v = a_dzs_deta_g
    dzs_dtau%v = a_dzs_dtau
    dzs_dxi%v = a_dzs_dxi
    dzs_dxi_g%v = a_dzs_dxi_g
    d_help_b%v = a_d_help_b
    d_help_c%v = a_d_help_c
    d_help_t%v = a_d_help_t
    ea = a_ea
    eaz_c = a_eaz_c
    eaz_c_quotient = a_eaz_c_quotient
#if (ENHMOD == 5)
    enh_c%v = a_enh_c
    enh_t%v = a_enh_t
#else
    enh_c = a_enh_c
    enh_t = a_enh_t
#endif
#if (CALCMOD >= 2)
    enth_c%v = a_enth_c
    enth_c_neu%v = a_enth_c_neu
#else
    enth_c = a_enth_c
    enth_c_neu = a_enth_c_neu
#endif
    enth_t = a_enth_t
    enth_t_neu = a_enth_t_neu
    et%v = a_et
    eta = a_eta
    F_INV = a_F_INV
    fc%v = a_fc
    flag_aa_nonzero = a_flag_aa_nonzero
    flag_calc_temp = a_flag_calc_temp
    flag_calving_front_1 = a_flag_calving_front_1
    flag_calving_front_2 = a_flag_calving_front_2
    flag_grounded_front_a_1 = a_flag_grounded_front_a_1
    flag_grounded_front_a_2 = a_flag_grounded_front_a_2
    flag_grounded_front_b_1 = a_flag_grounded_front_b_1
    flag_grounded_front_b_2 = a_flag_grounded_front_b_2
    flag_grounding_line_1 = a_flag_grounding_line_1
    flag_grounding_line_2 = a_flag_grounding_line_2
    flag_inner_point = a_flag_inner_point
    flag_shelfy_stream = a_flag_shelfy_stream
    flag_shelfy_stream_x = a_flag_shelfy_stream_x
    flag_shelfy_stream_y = a_flag_shelfy_stream_y
    flex_rig_lith = a_flex_rig_lith
#if (DYNAMICS==2 || MARGIN==3)
    flui_ave_sia%v = a_flui_ave_sia
#else
    flui_ave_sia = a_flui_ave_sia
#endif
    forcing_flag = a_forcing_flag
    G = a_G
    gamma_precip_lgm_anom = a_gamma_precip_lgm_anom
    gi_time_min = a_gi_time_min
    gi_time_max = a_gi_time_max
    gi_time_stp = a_gi_time_stp
#if TSURFACE==5
    allocate(glacial_index(0:a_ndata_gi))
    glacial_index = a_glacial_index
#endif
#if TSURFACE==4
    allocate(griptemp(0:a_ndata_grip))
    griptemp = a_griptemp
#endif
#if (DISC==2)
    glann_time_min = a_glann_time_min
    glann_time_stp = a_glann_time_stp
    glann_time_max = a_glann_time_max
    ndata_glann = a_ndata_glann
    dT_glann_CLIMBER = a_dT_glann_CLIMBER
#endif
    H_c%v = a_H_c
    H_c_neu%v = a_H_c_neu
    h_diff%v = a_h_diff
    H_R = a_H_R
    H_t%v = a_H_t
    H_target = a_H_target
    H_t_neu%v = a_H_t_neu
    H_w = a_H_w
    ij2n = a_ij2n
    insq_g11_g = a_insq_g11_g
    insq_g11_sgx = a_insq_g11_sgx
    insq_g22_g = a_insq_g22_g
    insq_g22_sgy = a_insq_g22_sgy
#if (defined(GRL) && DISC>0)
    iter_mar_coa_DW = a_iter_mar_coa_DW !ns
#endif
    KAPPA = a_KAPPA
    KAPPA_R = a_KAPPA_R
    KAPPA_C_IMQ = a_KAPPA_C_IMQ
    KAPPA_IMQ = a_KAPPA_IMQ !ns? called KAPPA, no _IMQ
    kc_cts = a_kc_cts
    kc_cts_neu = a_kc_cts_neu
    kei = a_kei
    kei_r_incr = a_kei_r_incr
    kei_r_max = a_kei_r_max
    L = a_L
    L_inv = a_L_inv !ns
    L_eto = a_L_eto !ns
    lambda = a_lambda
    LAMBDA0 = a_LAMBDA0
!#if OUTSER==3
    allocate(lambda_core(a_n_core))
    lambda_core = a_lambda_core ! there but is it working?
!#endif
#if (ENHMOD == 5)
    lambda_shear_c%v = a_lambda_shear_c
    lambda_shear_t%v = a_lambda_shear_t
#else
    lambda_shear_c = a_lambda_shear_c
    lambda_shear_t = a_lambda_shear_t
#endif
    maske = a_maske
    maske_maxextent = a_maske_maxextent
    maske_neu = a_maske_neu
    maske_old = a_maske_old
    maske_target = a_maske_target
    mask_ablation_type = a_mask_ablation_type
    n_slide_region = a_n_slide_region
    mb_source_apl = a_mb_source_apl
    accum_apl = a_accum_apl
    melt%v = a_melt
    melt_star%v = a_melt_star
    mult_test = a_mult_test
    MU_0 = a_MU_0
#if (defined(GRL) && DISC>0)
    m_D_DW = a_m_D_DW !ns
    m_H_DW = a_m_H_DW !ns
#endif
    n2i = a_n2i
    n2j = a_n2j
    ncid_core = a_ncid_core
    ncid_ser = a_ncid_ser
    ncid_temp_precip = a_ncid_temp_precip
    ndata_gi = a_ndata_gi
    ndata_grip = a_ndata_grip
    ndata_specmap = a_ndata_specmap
    ndata_temp_precip = a_ndata_temp_precip
    NUE = a_NUE
    n_core = a_n_core
    n_cts = a_n_cts
    n_cts_neu = a_n_cts_neu
#if (defined(GRL) && DISC>0)
    n_discharge_call_DW = a_n_discharge_call_DW !ns
#endif
    n_data_kei = a_n_data_kei !ns
    n_enth_min = a_n_enth_min !ns
    n_enth_max = a_n_enth_max !ns
#ifdef ANT
    n_bm_region = a_n_bm_region
#endif
    n_temp_min_IMQ = a_n_temp_min_IMQ !ns
    n_temp_max_IMQ = a_n_temp_max_IMQ !ns
    objf_test%v = a_objf_test
#if (CALCMOD >= 2)
    omega_c%v = a_omega_c
    omega_c_neu%v = a_omega_c_neu
#elif (CALCMOD == 0 || CALCMOD == 1)
    omega_c = a_omega_c
    omega_c_neu = a_omega_c_neu
#endif
    OMEGA_MAX = a_OMEGA_MAX
#if (CALCMOD >= 1) 
    omega_t%v = a_omega_t
    omega_t_neu%v = a_omega_t_neu
#else
    omega_t = a_omega_t
    omega_t_neu = a_omega_t_neu
#endif
    phi = a_phi
    PHI0 = a_PHI0
!#if OUTSER==3
    allocate(phi_core(a_n_core))
    phi_core = a_phi_core
!#endif
    PHI_SEP_0 = a_PHI_SEP_0
    PMAX_0 = a_PMAX_0
    precip_lgm_anom = a_precip_lgm_anom
    precip_ma_lgm_anom = a_precip_ma_lgm_anom
    precip_ma_present = a_precip_ma_present
    precip_present%v = a_precip_present
#if (defined(ANT) && SLIDE_LAW==2)
    p_b_w%v = a_p_b_w
#elif (SLIDE_LAW==1)
    p_b_w = a_p_b_w
#endif
    p_weert = a_p_weert
    p_weert_inv = a_p_weert_inv
#if (defined(ANT))
    qx%v = a_qx
    qy%v = a_qy
#else
    qx = a_qx
    qy = a_qy
#endif
    Q_bm%v = a_Q_bm
    Q_b_apl = a_Q_b_apl
    Q_b_tot%v = a_Q_b_tot
    q_geo%v = a_q_geo
    q_gl_g = a_q_gl_g
#if (CALCMOD >= 0) 
    Q_tld%v = a_Q_tld
#else
    Q_tld = a_Q_tld
#endif
    Q_w = a_Q_w
    Q_w_x = a_Q_w_x
    Q_w_y = a_Q_w_y
    q_weert = a_q_weert
    R = a_R
    rainfall%v = a_rainfall
#if (DYNAMICS==2 || MARGIN==3)
    ratio_sl_x%v = a_ratio_sl_x
    ratio_sl_y%v = a_ratio_sl_y
    ratio_sl%v   = a_ratio_sl
#else
    ratio_sl_x = a_ratio_sl_x
    ratio_sl_y = a_ratio_sl_y
    ratio_sl   = a_ratio_sl
#endif
    RF = a_RF
    RF_IMQ = a_RF_IMQ !just called RF in numcore
    RHO = a_RHO
    RHO_A = a_RHO_A
    RHO_C_IMQ = a_RHO_C_IMQ
    RHO_C_R = a_RHO_C_R
    RHO_I_IMQ = a_RHO_I_IMQ !ns
    RHO_SW = a_RHO_SW
    RHO_W = a_RHO_W
    runoff%v = a_runoff
    runoff_apl = a_runoff_apl
#if (defined(GRL) && DISC>0)
    r_mar_eff_DW = a_r_mar_eff_DW !ns
#endif
    R_T = a_R_T ! just called R_T
    R_T_IMQ = a_R_T_IMQ ! just called R_T
    sigma_c%v = a_sigma_c
    sigma_t%v = a_sigma_t
    smb_corr = a_smb_corr
    smb_corr_prescribed = a_smb_corr_prescribed    
    snowfall%v = a_snowfall
    specmap_time_max = a_specmap_time_max
    specmap_time_min = a_specmap_time_min
    specmap_time_stp = a_specmap_time_stp
#if (SEA_LEVEL==3)
    allocate(specmap_zsl(0:a_ndata_specmap))
    specmap_zsl = a_specmap_zsl !present but does header use this?
#endif
    sq_g11_g = a_sq_g11_g
    sq_g11_sgx = a_sq_g11_sgx
    sq_g11_sgy = a_sq_g11_sgy
    sq_g22_g = a_sq_g22_g
    sq_g22_sgx = a_sq_g22_sgx
    sq_g22_sgy = a_sq_g22_sgy
#if (defined(GRL) && DISC>0)
    s_dis_DW = a_s_dis_DW !ns
    c_dis_fac_DW = a_c_dis_fac_DW !ns
#endif
    S_STAT_0 = a_S_STAT_0
    target_topo_tau_0 = a_target_topo_tau_0
    temph_b = a_temph_b
    temp_b = a_temp_b
    temp_c%v = a_temp_c
    temp_c_m%v = a_temp_c_m
    temp_c_neu%v = a_temp_c_neu
    temp_ma_lgm_anom = a_temp_ma_lgm_anom
    temp_ma_present%v = a_temp_ma_present
    temp_mj_lgm_anom = a_temp_mj_lgm_anom
    temp_mj_present%v = a_temp_mj_present
    temp_mm_lgm_anom = a_temp_mm_lgm_anom
    temp_mm_present = a_temp_mm_present
#if ( (TSURFACE==6) && (ACCSURFACE==6) )
    allocate(temp_precip_time(0:a_ndata_temp_precip))
    temp_precip_time = a_temp_precip_time ! in there but does header file nee
#endif
    temp_precip_time_max = a_temp_precip_time_max
    temp_precip_time_min = a_temp_precip_time_min
    temp_precip_time_stp = a_temp_precip_time_stp
    temp_r%v = a_temp_r
    temp_r_neu%v = a_temp_r_neu
    temp_maat = a_temp_maat
    temp_s%v = a_temp_s
    temp_t_m%v = a_temp_t_m
    time_lag_asth = a_time_lag_asth
    time_target_topo_final = a_time_target_topo_final
    time_target_topo_init = a_time_target_topo_init
#if (defined(GRL) && DISC>0)
    T_sub_PD_DW = a_T_sub_PD_DW !ns
    T_sea_freeze_DW = a_T_sea_freeze_DW !ns
#endif
    txz_c = a_txz_c
    txz_t = a_txz_t
    tyz_c = a_tyz_c
    tyz_t = a_tyz_t
#if (DYNAMICS==2 || MARGIN==3)
    vis_ave_g%v = a_vis_ave_g
    vis_int_g%v = a_vis_int_g
#else
    vis_ave_g = a_vis_ave_g
    vis_int_g = a_vis_int_g
#endif
#if (DYNAMICS==2 || MARGIN==3)
    vx_b%v = a_vx_b
    vx_b_g%v = a_vx_b_g
#else
    vx_b = a_vx_b
    vx_b_g = a_vx_b_g
#endif
    vx_c%v = a_vx_c
    vx_g = a_vx_g
#if (defined(ANT))
    vx_m%v = a_vx_m
    vy_m%v = a_vy_m
    vx_m_sia%v = a_vx_m_sia
    vy_m_sia%v = a_vy_m_sia
    vx_m_ssa%v = a_vx_m_ssa
    vy_m_ssa%v = a_vy_m_ssa
#else
    vx_m = a_vx_m
    vy_m = a_vy_m
    vx_m_sia = a_vx_m_sia
    vy_m_sia = a_vy_m_sia
    vx_m_ssa = a_vx_m_ssa
    vy_m_ssa = a_vy_m_ssa
#endif
#if (DYNAMICS==2 || MARGIN==3)
    vx_s_g%v = a_vx_s_g
#else
    vx_s_g = a_vx_s_g
#endif
    vx_t%v = a_vx_t
#if (DYNAMICS==2 || MARGIN==3)
    vy_b%v = a_vy_b
    vy_b_g%v = a_vy_b_g
#else
    vy_b = a_vy_b
    vy_b_g = a_vy_b_g
#endif
    vy_c%v = a_vy_c
    vy_g = a_vy_g
#if (DYNAMICS==2 || MARGIN==3)
    vy_s_g%v = a_vy_s_g
#else
    vy_s_g = a_vy_s_g
#endif
    vy_t%v = a_vy_t
    vz_b%v = a_vz_b
    vz_c%v = a_vz_c
    vz_m%v = a_vz_m
#if (ENHMOD == 5 || DYNAMICS==2)
    vz_s%v = a_vz_s
#else
    vz_s = a_vz_s
#endif
    vz_t%v = a_vz_t
#if (REBOUND == 2)
    wss%v = a_wss
#else
    wss = a_wss
#endif
    xi = a_xi
#if OUTSER==3
    allocate(x_core(a_n_core))
    allocate(y_core(a_n_core))
    x_core = a_x_core
    y_core = a_y_core
#endif
    year_zero = a_year_zero
    zb%v = a_zb
    zb_neu%v = a_zb_neu
    zb_target = a_zb_target
    zeta_c = a_zeta_c
    zeta_r = a_zeta_r
    zeta_t = a_zeta_t
#if (REBOUND >= 1) 
    zl%v = a_zl
    zl_neu%v = a_zl_neu
#else
    zl = a_zl
    zl_neu = a_zl_neu
#endif
    zl0 = a_zl0
    zl_target = a_zl_target
    zm%v = a_zm
    zm_neu%v = a_zm_neu
    zs%v = a_zs
    zs_neu%v = a_zs_neu
    zs_ref = a_zs_ref
    zs_target = a_zs_target

  end subroutine var_transfer

!!-------------------------------------------------------------------------------
!!> This is another absolutely essential subroutine that is responsible for the
!!  writing of adjoint sensitivities to files for comparison against the gradient
!!  check as well as for other more science-targeted uses. As yet, it is rather
!!  ad-hoc and is structured for user-supplied controls. Every output file below is
!!  supplied and named directly by the user, and if for some reason the user runs
!!  an adjoint simulation where one of the controls is no longer ACTIVE, then
!!  errors may occur if syntaxes below are not proper. 
!!  For instance, if foobar is a control that is ACTIVE:
!!  in simulation (1) but not (2), then the proper syntax
!!  for the outputting of sensitivities (ignoring that outputting sensitivities for
!!  non-ACTIVE variables is a moot exercise) is:
!!   Simulation (1)
!!   "print stuff as below" foobar(j,i,k)%d
!!   Simulation (2)
!!   "print stuff as below" foobar(j,i,k) 
!!  Or, just comment out foobar from the writing to file entirely since all that
!!  will be printed is zeros anyway. And again, write to
!!  liz.curry.logan@gmail.com with error messages. 
!!<------------------------------------------------------------------------------
  subroutine print_output(runname) 

  use OAD_active
  use oad_sico_variables_m
  
  implicit none
  
  integer(i4b), parameter           :: points = 5
  integer(i4b), dimension(points)   :: ipoints, jpoints
  integer                           :: i, j, p
  character(len=100), intent(out)   :: runname
 
#if (!defined(ANT) && !defined(GRL))
  print *, "Adjoint mode only for Antarctica and Greenland right now."
#endif
 
  !-------- Open output files: 
  open(unit=79,file='AD_Vals_acc_fact.dat',status='replace')
  open(unit=80,file='AD_Vals_Q_tld.dat',status='replace')
  open(unit=81,file='AD_Vals_Q_bm.dat',status='replace')
  open(unit=82,file='AD_Vals_basal_stress.dat',status='replace')
  open(unit=83,file='AD_Vals_calving.dat',status='replace')
  open(unit=84,file='AD_Vals_dzs_deta_g.dat',status='replace')
  open(unit=85,file='AD_Vals_dzs_dxi_g.dat',status='replace')
  open(unit=86,file='AD_Vals_vis_int_g.dat',status='replace')
  open(unit=87,file='AD_Vals_precip_present_july.dat',status='replace')
  open(unit=88,file='AD_Vals_c_drag.dat',status='replace')
  open(unit=90,file='AD_Vals_q_geo.dat',status='replace')
  open(unit=93,file='AD_Vals_temp_c_base.dat',status='replace')
  open(unit=94,file='AD_Vals_temp_c_surf.dat',status='replace')
  open(unit=95,file='AD_Vals_H_c.dat',status='replace')
  open(unit=96,file='AD_Vals_for_grdchk.dat',status='replace')
  
  !-------- Outputting ad values for gradient-check-comparison:
  write(96, *) " (j,i)     control(j,i)"

  !-------- Select points along a spine of ice sheets
  do p = 1, points
#if (defined(GRL))
     ipoints(p) = int(real(IMAX/2))
     jpoints(p) = int(real(JMAX/5)) + (p-1) * points
#elif (defined(ANT))
     ipoints(p) = int(real(IMAX/3)) + int(real((.85-.33)*IMAX/points)) * (p - 1) 
     jpoints(p) = int(real(JMAX/2)) 
#endif
     write(96, '(f40.20)') H_c(jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') temp_c(KCMAX-3,jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') c_drag(jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') vx_c(24,jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') c_slide(jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') q_geo(jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') temp_ma_present(jpoints(p),ipoints(p))%d
     !write(96, '(f40.20)') precip_present(jpoints(p),ipoints(p),1)%d
     !write(96, '(f40.20)') acc_fact(jpoints(p),ipoints(p))%d
  end do
  

  !-------- Outputting ALL sensitivities to file (entire 2D fields):
  do j=0,JMAX
    do i=0,IMAX
  
      if ( isnan(H_c(j,i)%d) ) then
        write(unit=95,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=95,fmt='(t1,en18.6e4)',advance='no') H_c(j,i)%d
      end if
  
      if ( isnan(temp_c(KCMAX-3,j,i)%d) ) then
        write(unit=94,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=94,fmt='(t1,en18.6e4)',advance='no') temp_c(KCMAX-3,j,i)%d
      end if
  
      if ( isnan(temp_c(3,j,i)%d) ) then
        write(unit=93,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=93,fmt='(t1,en18.6e4)',advance='no') temp_c(3,j,i)%d
      end if
  
      if ( isnan(q_geo(j,i)%d) ) then
        write(unit=90,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=90,fmt='(t1,en18.6e4)',advance='no') q_geo(j,i)%d
      end if
  
#if (defined(ANT) && DYNAMICS==2 && MARGIN==3)
      if ( isnan(c_drag(j,i)%d) ) then
        write(unit=88,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=88,fmt='(t1,en18.6e4)',advance='no') c_drag(j,i)%d
      end if
#endif
  
      if ( isnan(precip_present(j,i,7)%d) ) then
        write(unit=87,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=87,fmt='(t1,en18.6e4)',advance='no') precip_present(j,i,7)%d
      end if
  
      if ( isnan(temp_ma_present(j,i)%d) ) then
        write(unit=93,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=93,fmt='(t1,en18.6e4)',advance='no') temp_ma_present(j,i)%d
      end if
  
  !------- Experimental scalar real for accumulation field
  !    if ( isnan(acc_fact(j,i)%d) ) then
  !      write(unit=79,fmt='(t1,en18.6e4)',advance='no') -66.6
  !    else
  !      write(unit=79,fmt='(t1,en18.6e4)',advance='no') acc_fact(j,i)%d
  !    end if
  
#if (defined(ANT) && DYNAMICS==2 && MARGIN==3)
      if ( isnan(temp_c(40,j,i)%d) ) then
        write(unit=86,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=86,fmt='(t1,en18.6e4)',advance='no') temp_c(40,j,i)%d
      end if
#endif 
  
      if ( isnan(dzs_dxi_g(j,i)%d) ) then
        write(unit=85,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=85,fmt='(t1,en18.6e4)',advance='no') dzs_dxi_g(j,i)%d
      end if
  
      if ( isnan(dzs_deta_g(j,i)%d) ) then
        write(unit=84,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=84,fmt='(t1,en18.6e4)',advance='no') dzs_deta_g(j,i)%d
      end if
  
#ifdef GRL
      if ( isnan(calving(j,i)%d) ) then
        write(unit=83,fmt='(t1,en18.6e4)',advance='no') -66.6
      else
        write(unit=83,fmt='(t1,en18.6e4)',advance='no') calving(j,i)%d
      end if
#endif
  
      if ( isnan(sigma_c(KCMAX,j,i)%d) ) then
        write(unit=82,fmt='(t1,en18.6e4)',advance='no') -66
      else
        write(unit=82,fmt='(t1,en18.6e4)',advance='no') sigma_c(KCMAX,j,i)%d
      end if
  
      if ( isnan(Q_bm(j,i)%d) ) then
        write(unit=81,fmt='(t1,en18.6e4)',advance='no') -66
      else
        write(unit=81,fmt='(t1,en18.6e4)',advance='no') Q_bm(j,i)%d
      end if
  
      if ( isnan(Q_tld(j,i)%d) ) then
        write(unit=80,fmt='(t1,en18.6e4)',advance='no') -66
      else
        write(unit=80,fmt='(t1,en18.6e4)',advance='no') Q_tld(j,i)%d
      end if
  
    end do
    write(unit=95,fmt='(t1)') 
    write(unit=94,fmt='(t1)') 
    write(unit=93,fmt='(t1)') 
    write(unit=90,fmt='(t1)') 
    write(unit=88,fmt='(t1)') 
    write(unit=87,fmt='(t1)') 
    write(unit=86,fmt='(t1)') 
    write(unit=85,fmt='(t1)') 
    write(unit=84,fmt='(t1)') 
    write(unit=83,fmt='(t1)') 
    write(unit=82,fmt='(t1)') 
    write(unit=81,fmt='(t1)') 
    write(unit=80,fmt='(t1)') 
  end do
  
  close(unit=96)
  close(unit=95)
  close(unit=94)
  close(unit=93)
  close(unit=90)
  close(unit=88)
  close(unit=87)
  close(unit=86)
  close(unit=85)
  close(unit=84)
  close(unit=83)
  close(unit=82)
  close(unit=81)
  close(unit=80)
  
  end subroutine print_output

#endif /* End ALLOW_OPENAD */

end module openad_m 
!
