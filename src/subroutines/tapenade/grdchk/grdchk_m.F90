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
    public  :: grdchk_main
  
  contains

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
   
   use cost_m
   
   implicit none
   
   integer(i4b)       :: ndat2d, ndat3d
   integer(i4b)       :: n_output
   real(dp)           :: delta_ts, glac_index
   real(dp)           :: mean_accum
   real(dp)           :: dtime, dtime_temp, dtime_wss, &
                                      dtime_out, dtime_ser
   real(dp)           :: time, time_init, time_end, time_output(100)
   real(dp)           :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
   real(dp)           :: z_mar
   
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
   open(99, file='GradientVals_'//trim(HEADER)//'.dat',&
       form="FORMATTED", status="REPLACE")
   open(98, file='CostVals_'//trim(HEADER)//'.dat',&
       form="FORMATTED", status="REPLACE")

!@ python_automated_grdchk IO begin @

	   open(9999,&
	   file='GradientVals_H_1.00E-03_repo_grl16_bm5_ss25ka_limited.dat',&
	   form="FORMATTED", status="REPLACE")
	   
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
        
            call sico_init(delta_ts, glac_index, &
                 mean_accum, &
                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 z_mar, &
                 ndat2d, ndat3d, n_output)

            perturbation = 1 + direction(d) * perturb_val 


          !-------- Controls to be perturbed (add your own here and below in
          !         subroutine print_output()
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir) 

                
            !@ python_automated_grdchk @

		            orig_val = H(j,i)
	                    H(j,i) = orig_val * perturbation
		
            ! Example -- H
            ! orig_val = H(j,i)
            ! H(j,i) = orig_val * perturbation

            ! -- sanity check
            write(6,fmt='(a,f40.20)') "orig_val = ", orig_val
            write(6,fmt='(a,f40.20)') "pert_val = ", orig_val*perturbation

            call sico_main_loop(delta_ts, glac_index, &
                 mean_accum, &
                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 z_mar, &
                 ndat2d, ndat3d, n_output)
          
            call cost_final()
            call sico_end

            ! Initialize compatible fields to 0
            ! 2D fields
            q_geo        = 0.0
            c_slide_init = 0.0
            H            = 0.0 ! Only compatible with ANF_DAT==1

            ! 3D fields
            temp_c       = 0.0 ! Not compatible with TEMP_INIT==5
       
            ! store cost
            fc_collected(d) = fc
        
            ! --------- calculate simple 2-sided finite difference due to
            !           perturbation: fc(+) - fc(-) / 2*espsilon
            if (orig_val .ne. 0) then
                gfd = (fc_collected(2) - fc_collected(3))/(2.d0 * perturb_val * orig_val)
            else
                gfd = 0.0
            end if          
          end do ! (close perturb loop)

          ! -- sanity check
          write(6, fmt='(a,f40.20)')   "Finite difference is = ", gfd
        
          ! --------- write these values to output file
          write(99, fmt='(f40.20)') gfd
          write(98, fmt='(f40.20)') fc_collected(1)
          write(98, fmt='(f40.20)') fc_collected(2)
          write(98, fmt='(f40.20)') fc_collected(3)
          write(98, fmt='(a)') '----------------------------------'
          
          !@ python_automated_grdchk IO write @
          write(9999, fmt='(f40.20)') gfd
    
   end do ! (close loop over points)
  
   close(unit=99)
   close(unit=98)
   !@ python_automated_grdchk IO end @
   close(unit=9999)
   end subroutine grdchk_main

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