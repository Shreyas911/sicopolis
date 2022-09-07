.. _AD:

Using Automatic Differentiation
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



Installation
============

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

Papers
======

.. _references:

.. only:: html

   .. rubric:: papers

Acknowledgements
================

This work was supported in part by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research Program, under contract number DE-AC02-06CH11357, and by National Science Foundation OPP/P2C2 grant #1903596.
