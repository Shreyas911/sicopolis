module ctrl_init_variables_m

  use sico_types_m
  use sico_variables_m
  use ctrl_map_genarr_m
  !use ctrl_map_gentim_m

  implicit none

  public :: ctrl_init_variables

contains

  subroutine ctrl_init_variables

    implicit none

    call ctrl_map_ini_genarr2d()
    !call ctrl_map_ini_genarr3d()
    !call ctrl_map_ini_gentim2d()

  end subroutine ctrl_init_variables

end module ctrl_init_variables_m
