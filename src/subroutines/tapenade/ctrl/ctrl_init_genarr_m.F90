#if defined(ALLOW_GENCTRL)
module ctrl_init_genarr_m

  use sico_types_m
  use sico_variables_m
  use ctrl_map_genarr_m

  implicit none

  public :: ctrl_init_genarr

contains

  subroutine ctrl_init_genarr

    implicit none

    call ctrl_map_ini_genarr2d()
    call ctrl_map_ini_genarr3d()

  end subroutine ctrl_init_genarr

end module ctrl_init_genarr_m
#endif
