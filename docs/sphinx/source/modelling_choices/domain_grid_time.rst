.. _domain_grid_time:

Model domain, grid and time
***************************

.. _model_domain:

Model domain
============

.. _defined_domain:

Selecting a pre-defined domain
------------------------------

SICOPOLIS provides several pre-defined model domains. They can be chosen by defining a domain code in the run-specs header as follows\:

.. code-block:: fortran

  #define GRL
  !                 Simulated domain:
  !                   ANT     - Antarctica
  !                   ASF     - Austfonna
  !                   EISMINT - EISMINT (Phase 2 SGE and modifications)
  !                   GRL     - Greenland
  !                   NHEM    - Northern hemisphere
  !                   SCAND   - Scandinavia
  !                   TIBET   - Tibet
  !                   NMARS   - North polar cap of Mars
  !                   SMARS   - South polar cap of Mars

This example would select the domain for the Greenland ice sheet. Correspondingly for the other listed domains.

.. _new_domain:

Setting up a new domain
-----------------------

In addition to the pre-defined domains, there is an unspecified domain XYZ. This framework allows creating new domains (Laurentide ice sheet, some ice cap, simple testing geometry...) quite easily. The directory src/subroutines/xyz, which hosts the domain-specific subroutines, is by default empty. If you want to create a new domain, copy the subroutines from the most similar existing domain, e.g., starting from Antarctica::

  cp src/subroutines/ant/*.F90 src/subroutines/xyz/

Then modify the routines according to your needs. Input files (topography etc.) must be placed in sico_in/xyz and specified in the run-specs header file \*.h as usual. The domain must be defined by the domain code

.. code-block:: fortran

  #define XYZ

in the header. For flexible testing, it is recommended to deactivate the compatibility check between horizontal resolution and number of grid points\:

.. code-block:: fortran

  #define CHECK_RES_IMAX_JMAX 0

If the new domain requires new global variables, they can be defined in the module src/subroutines/xyz/sico_vars.F90.

The subroutines for ISMIP HEINO (Calov et al. :cite:`calov_etal_2010`) are available in src/subroutines/xyz/heino, and the input files are in sico_in/xyz. If you copy the subroutines from src/subroutines/xyz/heino to src/subroutines/xyz, you can run ISMIP HEINO experiments (e.g., the run v5_heino50_st for which a run-specs header file is available).

.. _spatial_grid:

Spatial grid
============

Domain boundaries, resolution, projection...

Topography...

.. _model_time:

Model time
==========

Initial time, final time, time steps...
