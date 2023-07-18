.. _ad_running:

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

Here, ``{mode}`` refers to one of these options - ``normal, grdchk, adjoint, forward`` for normal (vanilla ``SICOPOLIS`` run), finite differences, adjoint and tangent linear mode respectively. ``{header}`` refers to the latter half of the name of the header file. If the header file is ``sico_specs_repo_grl20_ss25ka.h``, then the ``{header}`` is ``repo_grl20_ss25ka``. ``{domain}`` can either be ``grl`` or ``ant``, depending on Greenland or Antarctica. ``{dep_var}`` refers to the cost function / objective function. ``{ind_vars}`` is a list of independent variables for which we calculate the sensitivity of the objective function. ``{disc_and_grl}`` is either 0 or 1. It is 1 if the ``DISC`` flag in the header file is greater than 0 and the domain is Greenland.

The ``finite differences`` mode can only have one indpendent variable at a time, unlike the ``adjoint`` and ``forward`` modes.

The limitation here is that since Tapenade freshly generates the codes, the I/O for the differentiatied variables has to be done manually, at least for the current version of the setup. To tackle this, we have some Python utilities to automate this to a large extent.

Note that most of the Python utilities that have been developed for convenience have only been tested with one variable in ``{ind_vars}``. 
