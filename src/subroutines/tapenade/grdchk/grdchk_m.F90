!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  g r d c h k _ m
!
!> @file
!!
!! A catch-all module for grdchk-related subroutines. 
!!
!! @section Copyright
!!
!! Copyright 2017-2024 Shreyas Sunil Gaikwad,
!!                     Liz Curry-Logan, Sri Hari Krishna Narayanan,
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
!> Module for all tapenade-related subroutines 
!<------------------------------------------------------------------------------
#ifdef ALLOW_GRDCHK
module grdchk_m

    implicit none

    private :: deldirs
    public  :: grdchk_main, grdchk_array, grdchk_scalar
  
  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Subroutine :  g r d c h k _ m a i n
!   Purpose    :  Gradient check top level routine
!                 Compares the gradients calculated by the adjoint model
!                 to the gradients calculated by finite difference
!                 approximations
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine grdchk_main
  implicit none

!@ python_automated_grdchk scalar_or_array @
  
end subroutine grdchk_main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Subroutine :  g r d c h k _ s c a l a r
!   Purpose    :  Compares the gradients calculated by the adjoint model
!                 to the gradients calculated by finite difference
!                 approximations for 2D and 3D arrays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine grdchk_scalar
   
   use sico_types_m
   use sico_variables_m
#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
   use sico_vars_m
#endif
   
   use sico_init_m
   use sico_main_loop_m
   use sico_end_m
   
   use cost_m
   
   implicit none
   
   integer(i4b)       :: ndat2d, ndat3d
   integer(i4b)       :: n_output
   real(dp)           :: dtime, dtime_temp, dtime_wss, &
                                      dtime_out, dtime_ser
   real(dp)           :: time, time_init, time_end, time_output(100)
   real(dp)           :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
   
   !-------- Variable declarations needed for this routine specifically
   real(dp)                          :: orig_val, perturb_val = 0.001
   real(dp),     dimension(3)        :: fc_collected
   real(dp),     dimension(3)        :: direction
   real(dp)                          :: gfd0,gfd, perturbation
   integer(i4b)                      :: d
   character(len=100)                :: fname
   
   !-------- This array holds the direction of perturbations to follow:
   direction(1) = 0
   direction(2) = 1
   direction(3) = -1

#if (!defined(GRL) && !defined(ANT))
   print *, ">>> Adjoint only available for GRL and ANT right now; kill code." 
#endif

   !-------- Initialize output files 
   open(98, file='CostVals_'//trim(HEADER)//'.dat',&
       form="FORMATTED", status="REPLACE")

!@ python_automated_grdchk_scalar IO begin @

          !-------- Loop over perturbation direction (0, +, -)
          do d = 1, 3 

          !-------- Let yourself know where you are:
          print *, ' direction (d) [ ', d, ' ] '

          !-------- One complete forward run 
            call deldirs

            perturbation = 1 + direction(d) * perturb_val

          !-------- Controls to be perturbed before sico_init (add your own here)
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir)

            !@ python_automated_grdchk_scalar_before @

            ! Example -- gamma_s
            ! orig_val = gamma_s(1,1)
            ! if (orig_val .ne. 0) then
            !   gamma_s = orig_val * perturbation
            ! else
            !   gamma_s = perturbation-1
            ! end if

            call sico_init(dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 ndat2d, ndat3d, n_output)

          !-------- Controls to be perturbed after sico_init (add your own here)
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir)

            !@ python_automated_grdchk_scalar_after @

            ! Example -- gamma_s
            ! orig_val = gamma_s(1,1)
            ! if (orig_val .ne. 0) then
            !   gamma_s = orig_val * perturbation
            ! else
            !   gamma_s = perturbation-1
            ! end if

            ! -- sanity check
            write(6,fmt='(a,f40.20)') "orig_val = ", orig_val
            if (orig_val .ne. 0) then
              write(6,fmt='(a,f40.20)') "pert_val = ", orig_val * perturbation
            else
              write(6,fmt='(a,f40.20)') "pert_val = ", perturbation-1
            end if

            call sico_main_loop(dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 ndat2d, ndat3d, n_output)
          
            call cost_final()
            call sico_end

            ! Initialize compatible fields to 0
            ! 2D fields
            delta_tda_const = 0.0 ! Not compatible with TSURFACE>4
            q_geo           = 0.0
            c_slide_init    = 0.0
            p_weert         = 0.0
            q_weert         = 0.0
#if (ENHMOD==1 || ENHMOD==2 || ENHMOD==3)
            enh_fact_da_dummy2d_scalar = 0.0 ! Only used as a scalar ! Not compatible with ANF_DAT==3
#endif
#if (ENHMOD==2 || ENHMOD==3)
            enh_intg_da_dummy2d_scalar = 0.0 ! Only used as a scalar ! Not compatible with ANF_DAT==3
#endif
#if (FLOW_LAW==1)
            n_glen_da_dummy2d_scalar   = 0.0 ! Only used as a scalar ! Only compatible with :
                                             ! ((!defined(N_POWER_LAW_INT) OR N_POWER_LAW_INT < 0) && (defined(N_POWER_LAW_REAL)))
#endif
            H               = 0.0
            zs              = 0.0 ! Not compatible with ANF_DAT==2
            zl              = 0.0 ! Not compatible with ANF_DAT==2
            zl0             = 0.0
            zm              = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART) && CALCMOD!=1
            zb              = 0.0 ! Not compatible with ANF_DAT==2, or if ANF_DAT==1 and ZB_PRESENT_FILE is 'none' or not defined
#if (ACCSURFACE==2 || ACCSURFACE==3)
            gamma_s         = 0.0
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2)
            s_stat          = 0.0 ! Not compatible with SOLID_PRECIP==3, would need to also activate s_stat_solid_precip in that case.
            beta1           = 0.0
            beta2           = 0.0
            Pmax            = 0.0
            mu              = 0.0
#endif
#if (defined(GRL) && DISC>0)
            c_dis_da        = 0.0
#endif
#if (defined(PARAM_RHO_A))
            RHO_A           = 0.0
#endif
#if (REBOUND==1 || REBOUND==2)
            time_lag_asth   = 0.0
#endif
#if (REBOUND==2)
            flex_rig_lith   = 0.0
#endif
            ! 3D fields
            vx_c            = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART)
            vy_c            = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART)
            vz_c            = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART)
            temp_c          = 0.0 ! Not compatible with ANF_DAT==1 and TEMP_INIT==5
            age_c           = 0.0 ! Not compatible with ANF_DAT==1 and TEMP_INIT==5
            delta_tda       = 0.0 ! Not compatible with TSURFACE>4

            ! store cost
            fc_collected(d) = fc
        
          end do ! (close perturb loop)

          ! --------- calculate simple 2-sided finite difference due to
          !           perturbation: fc(+) - fc(-) / 2*espsilon
          if (orig_val .ne. 0) then
            gfd = (fc_collected(2) - fc_collected(3))/(2.d0 * perturb_val * orig_val)
          else
            gfd = (fc_collected(2) - fc_collected(1))/ perturb_val
          end if          

          ! -- sanity check
          write(6, fmt='(a,f40.20)')   "Finite difference is = ", gfd
        
          ! --------- write these values to output file
          write(99, fmt='(f40.20)') gfd
          write(98, fmt='(f40.20)') fc_collected(1)
          write(98, fmt='(f40.20)') fc_collected(2)
          write(98, fmt='(f40.20)') fc_collected(3)
          write(98, fmt='(a)') '----------------------------------'
          
          !@ python_automated_grdchk_scalar IO write @
  
   close(unit=98)
   !@ python_automated_grdchk_scalar IO end @
   end subroutine grdchk_scalar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Subroutine :  g r d c h k _ a r r a y
!   Purpose    :  Compares the gradients calculated by the adjoint model
!                 to the gradients calculated by finite difference
!                 approximations for 2D and 3D arrays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine grdchk_array
   
   use sico_types_m
   use sico_variables_m
#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
   use sico_vars_m
#endif
   use sico_init_m
   use sico_main_loop_m
   use sico_end_m
   
   use cost_m
   
   implicit none
   
   integer(i4b)       :: ndat2d, ndat3d
   integer(i4b)       :: n_output
   real(dp)           :: dtime, dtime_temp, dtime_wss, &
                                      dtime_out, dtime_ser
   real(dp)           :: time, time_init, time_end, time_output(100)
   real(dp)           :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
   
   !-------- Variable declarations needed for this routine specifically
   real(dp)                          :: orig_val, perturb_val = 0.001
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
   open(98, file='CostVals_'//trim(HEADER)//'.dat',&
       form="FORMATTED", status="REPLACE")

!@ python_automated_grdchk IO begin @
   
   !-------- Loop over points
   do p = 1, points !@ python_automated_grdchk limited_or_block_or_full @
     i = ipoints(p)
     j = jpoints(p)

          !-------- Loop over perturbation direction (0, +, -)
          do d = 1, 3 

          !-------- Let yourself know where you are:
          print *, ' point (p, i, j), direction (d) [ ', p , ', ', i, ', ', j, ', ', d, ' ] '

          !-------- One complete forward run 
            call deldirs

            perturbation = 1 + direction(d) * perturb_val

          !-------- Controls to be perturbed before sico_init (add your own here)
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir)

            !@ python_automated_grdchk_before @

            ! Example -- H
            ! orig_val = H(j,i)
            ! if (orig_val .ne. 0) then
            !   H(j,i) = orig_val * perturbation
            ! else
            !   H(j,i) = perturbation-1
            ! end if

            call sico_init(dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 ndat2d, ndat3d, n_output)

          !-------- Controls to be perturbed after sico_init (add your own here)
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir)

            !@ python_automated_grdchk_after @

            ! Example -- H
            ! orig_val = H(j,i)
            ! if (orig_val .ne. 0) then
            !   H(j,i) = orig_val * perturbation
            ! else
            !   H(j,i) = perturbation-1
            ! end if

            ! -- sanity check
            write(6,fmt='(a,f40.20)') "orig_val = ", orig_val
            if (orig_val .ne. 0) then
              write(6,fmt='(a,f40.20)') "pert_val = ", orig_val * perturbation
            else
              write(6,fmt='(a,f40.20)') "pert_val = ", perturbation-1
            end if

            call sico_main_loop(dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 ndat2d, ndat3d, n_output)
          
            call cost_final()
            call sico_end

            ! Initialize compatible fields to 0
            ! 2D fields
            delta_tda_const = 0.0 ! Not compatible with TSURFACE>4
            q_geo           = 0.0
            c_slide_init    = 0.0
            p_weert         = 0.0
            q_weert         = 0.0
#if (ENHMOD==1 || ENHMOD==2 || ENHMOD==3)
            enh_fact_da_dummy2d_scalar = 0.0 ! Only used as a scalar ! Not compatible with ANF_DAT==3
#endif
#if (ENHMOD==2 || ENHMOD==3)
            enh_intg_da_dummy2d_scalar = 0.0 ! Only used as a scalar ! Not compatible with ANF_DAT==3
#endif
#if (FLOW_LAW==1)
            n_glen_da_dummy2d_scalar   = 0.0 ! Only used as a scalar  ! Only compatible with :
                                             ! ((!defined(N_POWER_LAW_INT) OR N_POWER_LAW_INT < 0) && (defined(N_POWER_LAW_REAL)))
#endif
            H               = 0.0
            zs              = 0.0 ! Not compatible with ANF_DAT==2
            zl              = 0.0 ! Not compatible with ANF_DAT==2
            zl0             = 0.0
            zm              = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART) && CALCMOD!=1
            zb              = 0.0 ! Not compatible with ANF_DAT==2, or if ANF_DAT==1 and ZB_PRESENT_FILE is 'none' or not defined
#if (ACCSURFACE==2 || ACCSURFACE==3)
            gamma_s         = 0.0
#endif
#if (ABLSURFACE==1 || ABLSURFACE==2)
            s_stat          = 0.0 ! Not compatible with SOLID_PRECIP==3, would need to also activate s_stat_solid_precip in that case.
            beta1           = 0.0
            beta2           = 0.0
            Pmax            = 0.0
            mu              = 0.0
#endif
#if (defined(GRL) && DISC>0)
            c_dis_da        = 0.0
#endif
#if (defined(PARAM_RHO_A))
            RHO_A           = 0.0
#endif
#if (REBOUND==1 || REBOUND==2)
            time_lag_asth   = 0.0
#endif
#if (REBOUND==2)
            flex_rig_lith   = 0.0
#endif
            ! 3D fields
            vx_c            = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART)
            vy_c            = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART)
            vz_c            = 0.0 ! Only compatible with ANF_DAT==3 && !defined(LEGACY_RESTART)
            temp_c          = 0.0 ! Not compatible with ANF_DAT==1 and TEMP_INIT==5
            age_c           = 0.0 ! Not compatible with ANF_DAT==1 and TEMP_INIT==5
            delta_tda       = 0.0 ! Not compatible with TSURFACE>4
       
            ! store cost
            fc_collected(d) = fc

          end do ! (close perturb loop)

          ! --------- calculate simple 2-sided finite difference due to
          !           perturbation: fc(+) - fc(-) / 2*espsilon
          if (orig_val .ne. 0) then
            gfd = (fc_collected(2) - fc_collected(3))/(2.d0 * perturb_val * orig_val)
          else
            gfd = (fc_collected(2) - fc_collected(1))/ perturb_val
          end if          

          ! -- sanity check
          write(6, fmt='(a,f40.20)')   "Finite difference is = ", gfd
        
          ! --------- write these values to output file
          write(99, fmt='(f40.20)') gfd
          write(98, fmt='(f40.20)') fc_collected(1)
          write(98, fmt='(f40.20)') fc_collected(2)
          write(98, fmt='(f40.20)') fc_collected(3)
          write(98, fmt='(a)') '----------------------------------'
          
          !@ python_automated_grdchk IO write @
    
   end do ! (close loop over points)
  
   close(unit=98)
   !@ python_automated_grdchk IO end @
   end subroutine grdchk_array

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

end module grdchk_m
#endif /* ALLOW_GRDCHK */
