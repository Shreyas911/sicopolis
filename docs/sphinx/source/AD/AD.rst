.. _AD:

SICOPOLIS-AD v2
***************

Contributors
============

Shreyas Sunil Gaikwad, Laurent Hascoet, Sri Hari Krishna Narayanan, Liz Curry-Logan, Ralf Greve, Patrick Heimbach

Introduction and features
=========================

Previously, OpenAD have been used to get the adjoint of the SICOPOLIS code (`Logan et. al, 2020 <https://gmd.copernicus.org/articles/13/1845/2020/>`__). Our implementation has the following advantages over the previous implementation - 

1. It is up-to-date with the latest SICOPOLIS code

2. The AD tool Tapenade is open-source and actively maintained

3. A new tangent linear code generation capability is introduced (Forward Mode)

   * This is useful for adjoint validation, and from a physical perspective, also useful for Bayesian UQ and inverse modeling.
  
4. We are now able to deal with inputs in the NetCDF format

5. We have now correctly incorporated the external LIS solver, its tangent linear code, and its adjoint which improve the simulation of Antarctic ice shelves and Greenland outlet glaciers. (subject to more testing)

6. We leverage continuous integration and the pytest framework in order to track changes in the trunk that "break" the AD- based code generation - precludes the need for constant monitoring.

7. We “show” the entire code to Tapenade, including the initialization subroutines, thus avoiding cumbersome maintenance of subroutines OpenAD used to initialize for adjoint runs.

8. We have provided convenient Python scripts to make I/O with the differentiated variables easier.

Configuration
=============

.. toctree::
   :maxdepth: 2

   tapenade
   sicopolis_ad_config

Theory and tutorials
====================

.. toctree::
   :maxdepth: 2

   tutorials

SICOPOLIS-AD v2
===============

.. toctree::
   :maxdepth: 2

   running
   utilities

Acknowledgements
================
