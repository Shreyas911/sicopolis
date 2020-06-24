!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ e n d _ m
!
!> @file
!!
!! Ending of SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2020 Ralf Greve
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
!> Ending of SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_end_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of sico_end_m: Ending of SICOPOLIS.
!<------------------------------------------------------------------------------
  subroutine sico_end()

#if (NETCDF>1)
  use netcdf
  use nc_check_m
#endif

  implicit none

  integer(i4b) :: ierr

  character(len=64), parameter :: thisroutine = 'sico_end'

  write(unit=6, fmt='(/a/)') ' -------- sico_end --------'

  close(unit=12, status='keep')  ! Close time-series files
  close(unit=14, status='keep')
#if !defined(ALLOW_OPENAD)
  if (n_core >= 1) deallocate(lambda_core, phi_core, x_core, y_core, ch_core)
#endif

#if (defined(ASF) && WRITE_SER_FILE_STAKES>0)
  close(unit=41, status='keep')
  close(unit=42, status='keep')
  close(unit=43, status='keep')
  close(unit=44, status='keep')
  close(unit=45, status='keep')
  close(unit=46, status='keep')
  close(unit=47, status='keep')
  close(unit=48, status='keep')
  close(unit=49, status='keep')
  close(unit=50, status='keep')
#if !defined(ALLOW_OPENAD)
  deallocate(lambda_surf, phi_surf, x_surf, y_surf)
#endif
#endif

  if (allocated(specmap_zsl)) deallocate(specmap_zsl)

#if (defined(XYZ))
#if (defined(HEINO))
  close(unit=15, status='keep')
#endif
#endif

#if (NETCDF>1)
  call check( nf90_sync(ncid_ser),  thisroutine )
  call check( nf90_close(ncid_ser), thisroutine )
          ! Closing of NetCDF time-series output file
  if (n_core >= 1) then
     call check( nf90_sync(ncid_core),  thisroutine )
     call check( nf90_close(ncid_core), thisroutine )
  end if
          ! Closing of NetCDF time-series output file for the deep ice cores
#endif

#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
#if !defined(ALLOW_OPENAD)
  call lis_finalize(ierr)   ! Finalise execution environment of the
                            ! Library of Iterative Solvers Lis, if required
#else
  call lis_finalize_f(ierr) ! Finalise execution environment of the
                            ! Library of Iterative Solvers Lis, if required
#endif
#endif

  write(unit=6, fmt='(/a/)') '    * * * sicopolis.F90  r e a d y * * *'

  end subroutine sico_end

!-------------------------------------------------------------------------------

end module sico_end_m
!
