.. _AD:

Using automatic differentiation
*******************************

Contributors
============

* Shreyas Sunil Gaikwad - Principal developer and maintainer of SICOPOLIS-AD v2

* Laurent Hascoet - Developer of the open-source AD tool Tapenade

* Sri Hari Krishna Narayanan - Provided guidance on application of AD tools for both SICOPOLIS-AD v1 and SICOPOLIS-AD v2

* Liz Curry-Logan - Principal developer of SICOPOLIS-AD v1

* Ralf Greve - Developer of the SICOPOLIS ice sheet model

* Patrick Heimbach - PI of the NSF-supported project that funds this initiative

Introduction and features
=========================

Previously, OpenAD have been used to get the adjoint of the SICOPOLIS code (Logan et al. :cite:`logan_etal_2020`). The current implementation with Tapenade, SICOPOLIS-AD v2 (Gaikwad et al. :cite:`gaikwad_etal_2023`), has the following advantages over the previous implementation\:

1. It is up-to-date with the latest SICOPOLIS code

2. The AD tool Tapenade is open-source and actively maintained

3. A new tangent linear code generation capability is introduced (Forward Mode)

   * This is useful for adjoint validation, and from a physical perspective, also useful for Bayesian UQ and inverse modeling.
  
4. We are now able to deal with inputs in the NetCDF format

5. We have now correctly incorporated the external LIS solver, its tangent linear code, and its adjoint which improve the simulation of Antarctic ice shelves and Greenland outlet glaciers. (subject to more testing)

6. We leverage continuous integration and the pytest framework in order to track changes in the trunk that "break" the AD-based code generation - precludes the need for constant monitoring.

7. We "show" the entire code to Tapenade, including the initialization subroutines, thus avoiding cumbersome maintenance of subroutines OpenAD used to initialize for adjoint runs.

8. We have provided convenient Python scripts to make I/O with the differentiated variables easier.

In addition we also have the following previously available features\: 

1. The adjoint mode is available, like before along with the capability to do Finite Differences validation of the gradient computed using the adjoint mode.

The code has the following capabilities or possible applications\: 

1. Paleoclimatic inversions using the adjoint generated gradients as part of a model calibration exercise.

2. Uncertainty Quantification of calibrated parameters by leveraging the new tangent linear mode and adjoint mode.

3. Sensitivity analysis of various quantities of interest (QoI) to state parameters.

4. Optimal Experimental Design - where should sensors be placed such that the newly collected data optimally informs our uncertain parameters.

Installation
============

SICOPOLIS-AD v2 requires the installation of Tapenade as well as SICOPOLIS. It is mandatory to install the external libraries such as NetCDF, LIS to access the full functionality of the code, as well as git, to be able to clone and contribute to the repository.

.. toctree::
   :maxdepth: 2

   tapenade
   sicopolis_ad_config

.. _running:

Running SICOPOLIS-AD v2
=======================

The code for SICOPOLIS-AD v2 is kept mostly independent from the base SICOPOLIS code, allowing non-AD users to avoid it completely. All of the AD-related support routines and data files can be found in ``src/subroutines/tapenade``. Similarly, all utilities and testing files for AD simulations are stored in the ``test_ad`` directory. A separate Makefile is provided for AD purposes\: ``src/MakefileTapenade``.

``SICOPOLIS-AD v2`` can be run directly by interacting with the Makefile and the ``Fortran`` code. The user has to prepare a suitable header file for the base ``SICOPOLIS`` code and add a few more preprocessing directives to run the adjoint code with the same header file. This header file, along with a set of dependent and independent variables, is given to the Makefile. The Makefile executes the workflow for differentiating and compiling the code depending on the mode selected by the user (tangent linear, adjoint, finite differences, code coverage evaluation). The user must then insert the I/O statements in the Tapenade-generated code depending on what they wish to analyze. This is followed by recompilation and execution of the code.

Alternatively, the user can use the functions in ``test_ad/tapenade_config.py`` to automatically customize the setup. These functions are written for automated sensitivity studies. They can be easily modified to serve other purposes such as model calibration and UQ.

Typically, one runs the Makefile ``src/MakefileTapenade`` using the following generic command (assuming the current working directory is ``src/``)::

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
--------------------

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

Utilities
=========

.. toctree::
   :maxdepth: 2

   utilities

Theory and tutorials
====================

.. toctree::
   :maxdepth: 2

   theory
   tutorials

Acknowledgements
================

This work was supported in part by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research Program, under contract number DE-AC02-06CH11357, and by National Science Foundation OPP/P2C2 grant #1903596.
