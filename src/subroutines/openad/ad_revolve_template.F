#include "PACKAGES_CONFIG.h"

      subroutine template()
      use OAD_revolve_cp
      use OAD_cp
      use OAD_tape
      use OAD_rev
      use revolve

      
!$TEMPLATE_PRAGMA_DECLARATIONS
      LOGICAL :: initialized=.FALSE.
      TYPE(rvAction),save :: theAction
      CHARACTER(80) :: errorMsgRevolve
      integer, save :: jointCPCount
      integer, save :: currIter
      integer :: cp_loop_variable_1,cp_loop_variable_2,
     +     cp_loop_variable_3,cp_loop_variable_4,cp_loop_variable_5

      double precision, dimension(:), allocatable, save :: theArgFStack
      integer, save :: theArgFStackoffset=0, theArgFStackSize=0
      integer, dimension(:), allocatable, save :: theArgIStack
      integer, save :: theArgIStackoffset=0, theArgIStackSize=0
      type(modeType) :: our_orig_mode

      integer iaddr
      external iaddr
!      nIter0 = NINT( (startTime-baseTime)/deltaTClock )
      if (our_rev_mode%arg_store) then 
        call cp_write_open()
!$PLACEHOLDER_PRAGMA$ id=8
       call cp_close()
      end if 
      if (our_rev_mode%arg_restore) then
        call cp_read_open()
!$PLACEHOLDER_PRAGMA$ id=9
        call cp_close()
      end if
      if (our_rev_mode%plain) then
        do itercount=1, itercount_max
        CALL sico_main_loop_iter(delta_ts, glac_index,
     +                 mean_accum,
     +                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +                 time, time_init, time_end, time_output,
     +                 dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +                 z_sl, dzsl_dtau, z_mar,
     +                 ndat2d, ndat3d, n_output,
     +                 runname,
     +                 itercount,iter_temp,iter_wss,iter_ser,
     +                 iter_out,iter_output)
      enddo
      end if
      if (our_rev_mode%tape) then
         currIter=0
         jointCPcount=cp_fNumber()
         !initialized=rvInit(itercount_max,5,
         initialized=rvInit(itercount_max,120,
     +                      errorMsgRevolve,theAction)
         call rvVerbose(3)
         IF (.NOT.initialized) WRITE(*,'(A,A)') 'Error: ', errorMsgRevolve
         do while (theAction%actionFlag/=rvDone)
           theAction=rvNextAction()
           select case (theAction%actionFlag) 
             case (rvStore)
                call cp_write_open(theAction%cpNum+jointCPCount)
!$PLACEHOLDER_PRAGMA$ id=8
                call cp_close
             case (rvForward)
                call OAD_revPlain
                do currIter=currIter,theAction%iteration-1
                    CALL sico_main_loop_iter(delta_ts, glac_index,
     +               mean_accum,
     +               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +               time, time_init, time_end, time_output,
     +               dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +               z_sl, dzsl_dtau, z_mar,
     +               ndat2d, ndat3d, n_output,
     +               runname,
     +               currIter+1,iter_temp,iter_wss,iter_ser,
     +               iter_out,iter_output)
                end do
                call OAD_revTape
             case (rvFirstUTurn)
                      CALL sico_main_loop_iter(delta_ts, glac_index,
     +               mean_accum,
     +               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +               time, time_init, time_end, time_output,
     +               dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +               z_sl, dzsl_dtau, z_mar,
     +               ndat2d, ndat3d, n_output,
     +               runname,
     +               currIter+1,iter_temp,iter_wss,iter_ser,
     +               iter_out,iter_output)
! get out now ...
                exit                
             end select   
          end do
      end if 
      if (our_rev_mode%adjoint) then
        IF (.NOT.initialized) WRITE(*,'(A)') 'Error: not initialized'
        do while (theAction%actionFlag/=rvDone)
           select case (theAction%actionFlag) 
             case (rvFirstUTurn)
!we taped already ... see above
                      CALL sico_main_loop_iter(delta_ts, glac_index,
     +               mean_accum,
     +               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +               time, time_init, time_end, time_output,
     +               dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +               z_sl, dzsl_dtau, z_mar,
     +               ndat2d, ndat3d, n_output,
     +               runname,
     +               currIter+1,iter_temp,iter_wss,iter_ser,
     +               iter_out,iter_output)
             case (rvStore)
                call cp_write_open(theAction%cpNum+jointCPCount)
!$PLACEHOLDER_PRAGMA$ id=8
                call cp_close
             case (rvRestore)
                call cp_read_open(theAction%cpNum+jointCPCount)
!$PLACEHOLDER_PRAGMA$ id=9
                currIter=theAction%iteration
                call cp_close
             case (rvForward)
                call OAD_revPlain
                do currIter=currIter,theAction%iteration-1
                         CALL sico_main_loop_iter(delta_ts, glac_index,
     +               mean_accum,
     +               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +               time, time_init, time_end, time_output,
     +               dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +               z_sl, dzsl_dtau, z_mar,
     +               ndat2d, ndat3d, n_output,
     +               runname,
     +               currIter+1,iter_temp,iter_wss,iter_ser,
     +               iter_out,iter_output)
                end do
                call OAD_revAdjoint
             case (rvUTurn)
                call OAD_revTape
                      CALL sico_main_loop_iter(delta_ts, glac_index,
     +               mean_accum,
     +               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +               time, time_init, time_end, time_output,
     +               dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +               z_sl, dzsl_dtau, z_mar,
     +               ndat2d, ndat3d, n_output,
     +               runname,
     +               itercount,iter_temp,iter_wss,iter_ser,
     +               iter_out,iter_output)
                call OAD_revAdjoint
                      CALL sico_main_loop_iter(delta_ts, glac_index,
     +               mean_accum,
     +               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser,
     +               time, time_init, time_end, time_output,
     +               dxi, deta, dzeta_c, dzeta_t, dzeta_r,
     +               z_sl, dzsl_dtau, z_mar,
     +               ndat2d, ndat3d, n_output,
     +               runname,
     +               currIter+1,iter_temp,iter_wss,iter_ser,
     +               iter_out,iter_output)
           end select   
           theAction=rvNextAction()
        end do
      end if 
      end subroutine template

