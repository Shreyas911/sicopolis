.. _running:

Running SICOPOLIS-AD v2
***********************

SICOPOLIS-AD v2 runs almost independently of the setup for SICOPOLIS. It has its own Makefile ``src/MakefileTapenade``, typically one runs it using the following generic command (assuming the current working directory is ``src/``)::

    % make -f MakefileTapenade clean
    % make -f MakefileTapenade driver{mode} HEADER={header} DOMAIN_SHORT={domain} DEP_VAR={dep_var} IND_VARS={ind_vars} DISC_AND_GRL={disc_and_grl}
    #### Add I/O commands in the differentiated code if needed and recompile
    % make -f MakefileTapenade driver{mode} HEADER={header} DOMAIN_SHORT={domain} DEP_VAR={dep_var} IND_VARS={ind_vars} DISC_AND_GRL={disc_and_grl}
    % ./driver{mode}

Here, ``{mode}`` refers to one of these options - ``normal, grdchk, adjoint, forward`` for normal (vanilla ``SICOPOLIS`` run), finite differences, adjoint and tangent linear mode respectively. ``{header}`` refers to the latter half of the name of the header file. If the header file is ``sico_specs_v5_grl20_ss25ka.h``, then the ``{header}`` is ``v5_grl20_ss25ka``. ``{domain}`` can either be ``grl`` or ``ant``, depending on Greenland or Antarctica. ``{dep_var}`` refers to the cost function / objective function. ``{ind_vars}`` is a list of independent variables for which we calculate the sensitivity of the objective function. ``{disc_and_grl}`` is either 0 or 1. It is 1 if the ``DISC`` flag in the header file is greater than 0 and the domain is Greenland.

The ``finite differences`` mode can only have one indpendent variable at a time, unlike the ``adjoint`` and ``forward`` modes.

The limitation here is that since Tapenade freshly generates the codes, the I/O for the differentiatied variables has to be done manually, at least for the current version of the setup. To tackle this, we have some Python utilities to automate this to a large extent.


Note that most of the Python utilities that have been developed for convenience have only been tested with one variable in ``{ind_vars}``. 

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


