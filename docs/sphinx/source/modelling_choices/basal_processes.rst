.. _subglacial_processes:

Subglacial processes
********************

.. _basal_sliding:

Basal sliding
=============

Lorem ipsum...

.. _ghf:

Geothermal heat flux
====================

The geothermal heat flux, assumed to be time-independent, can be specified in the run-specs header as either a constant value or a spatially varying distribution. For the example of simulation :file:`v5_grl16_bm5_ss25ka`:

.. code-block:: fortran

  #define Q_GEO_MOD 2
  !                       1 : Constant geothermal heat flux defined
  !                           by parameter Q_GEO
  !                       2 : Spatially varying geothermal heat flux
  !                           read from file

  #define Q_GEO 0.0d0
  !                       Constant geothermal heat flux (for Q_GEO_MOD==1),
  !                       in mW/m2

  #define Q_GEO_FILE 'GHF_Greenland_Ver2.0_GridEPSG3413_16km.nc'
  !                       Name of the file containing the spatially varying
  !                       geothermal heat flux (for Q_GEO_MOD==2)
  !                       (if NetCDF, variable name must be 'GHF')

If, like in the above example, a file with gridded data is provided, it must match the chosen horizontal grid (see Section ":ref:`spatial_grid`"). The format can either be NetCDF (:file:`*.nc`) or ASCII (any other file extension).

A further, relevant parameter is ``Q_LITHO`` from the "Lithosphere (bedrock) modelling" section in the run-specs header:

.. code-block:: fortran

  #define Q_LITHO 0
  !                       0 : No coupled heat-conducting bedrock
  !                       1 : Coupled heat-conducting bedrock

If set to ``0``, the geothermal heat flux is imposed directly at the grounded ice base, which is suitable for steady-state simulations because it reduces the time required to reach the steady state. However, for transient simulations, ``1`` is the preferred setting. The geothermal heat flux is then imposed at the base of the lithosphere layer (thickness defined by variable ``H_R`` in the physical-parameter file, see Section ":ref:`getting_started-phys_para`"), so that the thermal inertia of the lithosphere is properly accounted for.

.. _basal_hydrology:

Basal hydrology
===============

Basal hydrology can be selected in the run-specs header by the parameter ``BASAL_HYDROLOGY``\:

* If set to ``0``, basal hydrology is ignored.

* If set to ``1``, it is assumed that basal water exists and moves in a thin (order of millimetres) and distributed water film. The film thickness is computed by a steady-state routing scheme for subglacial water that receives its input from the basal melting rate under grounded ice (Le Brocq et al. :cite:`lebrocq_etal_2006, lebrocq_etal_2009`, Calov et al. :cite:`calov_etal_2018`). The computations are carried out by the module ``hydro_m``.

.. _gia:

Glacial isostatic adjustment
============================

Lorem ipsum...
