.. _clim_ocn_forcing:

Oceanic forcing
***************

.. _sea_level:

Sea level
=========

The sea level surrounding the simulated ice sheet determines the land area available for glaciation. The level :math:`z=0` corresponds to the mean sea level at present day, and the SICOPOLIS variable ``z_sl`` denotes the sea level relative to this reference.

Two different options for prescribing the sea level are available, selected in the run-specs headers by the parameter ``SEA_LEVEL``\:

* ``1``: Temporally constant sea level z_sl, specified by the parameter ``Z_SL0``.

* ``3``: Time-dependent sea level (e.g., reconstruction from data) from an input file (ASCII or NetCDF), specified by the parameter ``SEA_LEVEL_FILE``.

Both options assume a spatially constant sea level. However, the variable ``z_sl`` is actually a 2D field, so that SICOPOLIS can handle in principle a spatially variable sea level as well.

.. _ice_shelf_basal_melting:

Ice-shelf basal melting
=======================

The parameter ``FLOATING_ICE_BASAL_MELTING`` in the run-specs headers allows specifying the melting rate under ice shelves (floating ice). For all terrestrial ice sheets, the following options can be chosen\:

* ``1``: Constant values for the continental shelf and the abyssal ocean, respectively.

* ``4``: Local parameterization as a function of the oceanic thermal forcing (ocean temperature :math:`T_\mathrm{oc}` minus ice shelf basal temperature :math:`T_\mathrm{m,b}`).

The parameterization for ``FLOATING_ICE_BASAL_MELTING = 4`` reads

.. math::
  :label: ice_shelf_bas_melt_1

  a_\mathrm{b} = S_\mathrm{w} \Omega\,(T_\mathrm{oc}-T_\mathrm{m,b})^\alpha\,,

where :math:`a_\mathrm{b}` is the sub-ice-shelf melting rate. The parameters :math:`\Omega`, :math:`T_\mathrm{oc}` and :math:`\alpha` can be set in the run-specs header. If the thermal forcing is negative, it is set to zero. The scaling factor :math:`S_\mathrm{w}` is

.. math::
  :label: ice_shelf_bas_melt_scaling_factor

  S_\mathrm{w}
    = \mathrm{tanh}\,\bigg(\frac{H_\mathrm{w}}{H_\mathrm{w,0}}\bigg)\,.

This factor reduces the melting rate close to the grounding line where the water column :math:`H_\mathrm{w}` is thin. The parameter :math:`H_\mathrm{w,0}` can also be set in the run-specs header; a value recommended by Asay-Davis et al. :cite:`asay-davis_etal_2016` is :math:`75\,\mathrm{m}`. Setting this parameter to zero results in :math:`S_\mathrm{w}=1` everywhere; the scaling is then switched off.

For the Antarctic ice sheet, two additional options are available\:

* ``5``: Sector-wise, local parameterization as a function of the thermal forcing (Greve and Galton-Fenzi :cite:`greve_galton-fenzi_2017`).

* ``6``: "ISMIP6 standard approach": Sector-wise, non-local quadratic parameterization as a function of the thermal forcing (Jourdain et al. :cite:`jourdain_etal_2020`).

Currently, option ``FLOATING_ICE_BASAL_MELTING = 6`` is the only one that allows prescribing a time-dependent thermal forcing.

  .. _calving_ice_shelves:

Ice-shelf calving
=================

The options for calving of ice shelves (floating ice) can be selected in the run-specs headers by the parameter ``ICE_SHELF_CALVING``\:

* ``1``: Unlimited expansion of ice shelves, no calving.

* ``2``: Instantaneous calving of ice shelves if the thickness is less than a threshold thickness, specified by the parameter ``H_CALV``.

* ``3``: "Float-kill": Instantaneous removal of all floating ice.

For the Antarctic ice sheet, yearly ISMIP6-type ice-shelf collapse masks can be prescribed (Seroussi et al. :cite:`seroussi_etal_2020`). This requires the setting ``ICE_SHELF_COLLAPSE_MASK = 1`` and additional parameters as described in the run-specs headers.

.. _calving_marine_ice:

Marine-ice calving
==================

For calving of grounded marine ice, the following options are available\:

* | Parameterization for "underwater-ice" calving (Dunse et al. :cite:`dunse_etal_2011`).
  | To be selected by the following combination of run-specs-header parameters\:
  | ``MARGIN = 2``, ``MARINE_ICE_FORMATION = 2``, ``MARINE_ICE_CALVING = 9``.
  | Further parameters (``CALV_UW_COEFF``, ``R1_CALV_UW``, ``R2_CALV_UW``) as described in the run-specs headers.

For the Greenland ice sheet, yearly ISMIP6-type retreat masks can be prescribed (Goelzer et al. :cite:`goelzer_etal_2020`). This requires the setting ``RETREAT_MASK = 1`` and additional parameters as described in the run-specs headers.
