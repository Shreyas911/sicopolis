.. _clim_ocn_forcing:

Oceanic forcing
***************

.. |nbsp| unicode:: 0xA0 
   :trim:

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

* ``1``: Constant values for the continental shelf (parameter ``QBM_FLOAT_1``) and the abyssal ocean (``QBM_FLOAT_3``), respectively. The threshold seabed elevation separating the continental shelf from the abyssal ocean is defined by the parameter ``Z_ABYSS``.

* ``4``: Local parameterization as a function of the oceanic thermal forcing :math:`T_\mathrm{f}=T_\mathrm{oc}-T_\mathrm{b}` (difference between the ocean temperature :math:`T_\mathrm{oc}` and the ice-shelf basal temperature :math:`T_\mathrm{b}`)\:

  .. math::
    :label: ice_shelf_bas_melt_1

    a_\mathrm{b} = \Omega\,T_\mathrm{f}^\alpha\,,

  where :math:`a_\mathrm{b}` is the sub-ice-shelf melting rate. The parameters :math:`T_\mathrm{oc}`, :math:`\Omega` and :math:`\alpha` can be set in the run-specs header (``TEMP_OCEAN``, ``OMEGA_QBM`` and ``ALPHA_QBM``, respectively). The ice-shelf basal temperature is computed as

  .. math::
    :label: ice_shelf_bas_temp

    T_\mathrm{b} = -\beta_\mathrm{sw} d 
                     - \Delta{}T_\mathrm{m,sw}\,,

  where :math:`d` is the draft (depth of the ice-shelf base below sea level), :math:`\beta_\mathrm{sw}=7.61\times{}10^{-4}\,\mathrm{K\,m^{-1}}` the Clausius-Clapeyron gradient and :math:`\Delta{}T_\mathrm{m,sw}=1.85^\circ\mathrm{C}` the melting-point lowering due to the average salinity of sea water.

For the Antarctic ice sheet, two additional options are available\:

* ``5``: Sector-wise, local parameterization as a function of the thermal forcing (Greve and Galton-Fenzi :cite:`greve_galton-fenzi_2017`). This parameterization is modified after Beckmann and Goosse :cite:`beckmann_goosse_2003`, with a linear dependence on the thermal forcing :math:`T_\mathrm{f}` and an additional power-law dependence on the draft :math:`d`\:

  .. math::
    :label: ice_shelf_bas_melt_2

    a_\mathrm{b} 
      = \frac{\rho_\mathrm{sw}c_\mathrm{sw}\gamma_\mathrm{t}}{\rho L}
        \,\Omega\,\bigg(\frac{d}{d_0}\bigg)^\alpha \,T_\mathrm{f}\,,

  where :math:`\rho` and :math:`\rho_\mathrm{sw}` are the ice and sea-water densities, :math:`L` is the latent heat of melting (all defined in the physical-parameter file), :math:`c_\mathrm{sw}=3974\,\mathrm{J\,kg^{-1}\,K^{-1}}` is the specific heat of sea water, :math:`\gamma_\mathrm{t}=5\times{}10^{-5}\,\mathrm{m\,s^{-1}}` is the exchange velocity for temperature and :math:`d_0=200\,\mathrm{m}` is the reference draft.

  The parameters :math:`\Omega` and :math:`\alpha` result from a tuning procedure for eight different sectors, using observed present-day melt rates as a target (as explained in the main part and appendix of Greve and Galton-Fenzi :cite:`greve_galton-fenzi_2017`). For the thermal forcing :math:`T_\mathrm{f}`, :math:`T_\mathrm{oc}` is chosen for each sector as the sector-averaged temperature at 500 metres depth just outside the ice-shelf cavity (computed with data from the World Ocean Atlas 2009), while :math:`T_\mathrm{b}` is computed by Eq. |nbsp| :eq:`ice_shelf_bas_temp`.

* ``6``: "ISMIP6 standard approach": Sector-wise, non-local quadratic parameterization as a function of the thermal forcing (Jourdain et al. :cite:`jourdain_etal_2020`, Seroussi et al. :cite:`seroussi_etal_2020`).

  Currently, this is the only option that allows prescribing a time-dependent thermal forcing. Following the ISMIP6-Antarctica protocol, it must be provided as NetCDF input files that contain for each year the mean-annual, 3D thermal forcing for the entire computational domain.

For all cases, an additional scaling factor :math:`S_\mathrm{w}` can be applied (:math:`a_\mathrm{b}\rightarrow{}S_\mathrm{w}\,a_\mathrm{b}`), defined as

.. math::
  :label: ice_shelf_bas_melt_scaling_factor

  S_\mathrm{w}
    = \mathrm{tanh}\,\bigg(\frac{H_\mathrm{w}}{H_\mathrm{w,0}}\bigg)\,.

This factor reduces the melting rate close to the grounding line where the water column :math:`H_\mathrm{w}` is thin. The parameter :math:`H_\mathrm{w,0}` can be set in the run-specs header (``H_W_0``). A value recommended by Asay-Davis et al. :cite:`asay-davis_etal_2016` is :math:`75\,\mathrm{m}`, while Gladstone et al. :cite:`gladstone_etal_2017` used :math:`36.79\,(=100/e)\,\mathrm{m}`. Setting this parameter to zero results in :math:`S_\mathrm{w}=1` everywhere; the scaling is then switched off.

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
