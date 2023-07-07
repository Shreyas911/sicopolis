.. _ad_activating_variables:

Activating variables
====================

There are some variables in the SICOPOLIS code, for example ``c_slide``, the basal sliding coefficient, that are constant in time and reinitialized at each time step in the subroutines in ``src/subroutines/general/calc_vxy_m.F90``. This reinitialization has an "erasing" effect on the the dependency of cost functions to such variables, resulting in identically zero adjoint fields. From a physics standpoint, it is pretty clear that ``c_slide`` should have affect various kinds of cost functions.

The workaround to this issue is to complete the following steps (the code snippets are based on the ``c_slide`` example)\:

1. Declare a dummy variable in ``src/subroutines/general/sico_variables_m.F90``.

  .. code-block:: fortran

    real(dp), dimension(0:JMAX,0:IMAX) :: xxc_slide

2. Initialize the dummy variable to 0. at the end of ``src/subroutines/{domain}/sico_init_m.F90``. Here ``{domain}`` can either be ``grl`` or ``ant``.

  .. code-block:: fortran

    xxc_slide = 0.

3. Add the dummy variable to the actual variable. In the case of ``c_slide``, add it at the point where it's initialization is complete. 

  .. code-block:: fortran

    c_slide = c_slide + xxc_slide

4. Evaluate sensitivities to the dummy variable i.e. the ``{ind_vars}`` should contain ``xxc_slide`` for a non-trivial senstivity with respect to ``c_slide``.

This procedure has the effect of continuity of the variable ``xxc_slide`` in time, preventing the "erasing" effect caused by the reinitialization of ``c_slide`` to the same value for each time step. Also, since ``xxc_slide`` is initialized to 0, the results in the normal simulations are unaffected, because for the purposes of the non-linear forward model, all we did was add zeros to ``c_slide``.
