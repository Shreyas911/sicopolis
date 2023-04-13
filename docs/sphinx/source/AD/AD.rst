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

How to contribute?
==================

SICOPOLIS-AD v2 is an open source project that relies on the participation of its users, and we welcome contributions. Users can contribute using the usual pull request mechanisms in git, and if the contribution is substantial, they can contact us to discuss gaining direct access to the repository.

If you think you’ve found a bug, please check if you’re using the latest version of the model. If the bug is still present, then think about how you might fix it and file a ticket in the Gitlab issue tracker (you might need to request membership access on Gitlab, which we can approve). Your ticket should include: what the bug does, the location of the bug: file name and line number(s), and any suggestions you have for how it might be fixed.

To request a new feature, or guidance on how to implement it yourself, please open a ticket with a clear explanation of what the feature will do.

You can also directly contact Shreyas Gaikwad (shreyas.gaikwad@utexas.edu) for any of the above.

Introduction and features
=========================

Previously, OpenAD have been used to get the adjoint of the SICOPOLIS code (`Logan et. al, 2020 <https://gmd.copernicus.org/articles/13/1845/2020/>`__). The current implementation with Tapenade (SICOPOLIS-AD v2) has the following advantages over the previous implementation - 

1. It is up-to-date with the latest SICOPOLIS code

2. The AD tool Tapenade is open-source and actively maintained

3. A new tangent linear code generation capability is introduced (Forward Mode)

   * This is useful for adjoint validation, and from a physical perspective, also useful for Bayesian UQ and inverse modeling.
  
4. We are now able to deal with inputs in the NetCDF format

5. We have now correctly incorporated the external LIS solver, its tangent linear code, and its adjoint which improve the simulation of Antarctic ice shelves and Greenland outlet glaciers. (subject to more testing)

6. We leverage continuous integration and the pytest framework in order to track changes in the trunk that "break" the AD- based code generation - precludes the need for constant monitoring.

7. We “show” the entire code to Tapenade, including the initialization subroutines, thus avoiding cumbersome maintenance of subroutines OpenAD used to initialize for adjoint runs.

8. We have provided convenient Python scripts to make I/O with the differentiated variables easier.

In addition we also have the following previously available features - 

1. The adjoint mode is available, like before along with the capability to do Finite Differences validation of the gradient computed using the adjoint mode.

The code has the following capabilities or possible applications - 

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

Theory and tutorials
====================

.. toctree::
   :maxdepth: 2

   theory
   tutorials

SICOPOLIS-AD v2
===============

The code for SICOPOLIS-AD v2 is kept mostly independent from the base SICOPOLIS code, allowing non-AD users to avoid it completely. All of the AD-related support routines and data files can be found in ``src/subroutines/tapenade``. Similarly, all utilities and testing files for AD simulations are stored in the ``test_ad`` directory. A separate Makefile is provided for AD purposes - ``src/MakefileTapenade``.

``SICOPOLIS-AD v2`` can be run directly by interacting with the Makefile and the ``Fortran`` code. The user has to prepare a suitable header file for the base ``SICOPOLIS`` code and add a few more preprocessing directives to run the adjoint code with the same header file. This header file, along with a set of dependent and independent variables, is given to the Makefile. The Makefile executes the workflow for differentiating and compiling the code depending on the mode selected by the user (tangent linear, adjoint, finite differences, code coverage evaluation). The user must then insert the I/O statements in the Tapenade-generated code depending on what they wish to analyze. This is followed by recompilation and execution of the code.

Alternatively, the user can use the functions in ``test_ad/tapenade_config.py`` to automatically customize the setup. These functions are written for automated sensitivity studies. They can be easily modified to serve other purposes such as model calibration and UQ.

.. toctree::
   :maxdepth: 2

   running
   utilities

Acknowledgements
================

This work was supported in part by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research Program, under contract number DE-AC02-06CH11357, and by National Science Foundation OPP/P2C2 grant #1903596.
