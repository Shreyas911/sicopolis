module ctrl_init_gen_m

  use sico_types_m
  use sico_variables_m
  use ctrl_map_genarr_m
  use ctrl_map_gentim_m

  implicit none

  public :: ctrl_init_gen

contains

  subroutine ctrl_init_gen

    implicit none

#ifdef DO_CTRL_GENARR2D
    call ctrl_map_ini_genarr2d()
#endif
#ifdef DO_CTRL_GENARR3D
    call ctrl_map_ini_genarr3d()
#endif
#ifdef DO_CTRL_GENTIM2D
    call ctrl_map_ini_gentim2d()
#endif

  end subroutine ctrl_init_gen

end module ctrl_init_gen_m
