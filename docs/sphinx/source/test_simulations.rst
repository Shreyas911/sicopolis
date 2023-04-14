.. _test_simulations:

Test simulations
****************

A number of computationally rather inexpensive test runs...

     * ``v5_vialov3d25``

       * 3-d version of the 2-d Vialov profile

       * SIA, resolution 25 km, :math:`t=0\ldots{}100\,\mathrm{ka}`.

       * Similar to the EISMINT Phase 1 fixed-margin experiment

       * Huybrechts et al. :cite:`huybrechts_etal_1996`, but without thermodynamics. Instead, isothermal conditions with :math:`T=-10^{\circ}` C everywhere are assumed.

     * ``v5_emtp2sge25_expA``

       * EISMINT Phase 2 Simplified Geometry Experiment A

       * SIA, resolution 25 km, :math:`t=0\ldots{}200\,\mathrm{ka}` (Payne et al. :cite:`payne_etal_2000`).

       * The thermodynamics solver for this run is the one-layer melting-CTS enthalpy scheme (ENTM), while all other runs employ the polythermal two-layer scheme (POLY) (Greve and Blatter :cite:`greve_blatter_2016`).

     * ``v5_grl16_bm5_ss25ka``

       * Greenland ice sheet, SIA, resolution 16 km, short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for modern climate conditions (unpublished).

     * ``v5_ant40_b2_ss25ka`` 

       * Antarctic ice sheet without ice shelves, SIA, resolution 40 km, short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for modern climate conditions (unpublished).

     * ``v5_grl20_b2_paleo21``

       * Greenland ice sheet, SIA, resolution 20 km, :math:`t=-140\ldots{}0\,\mathrm{ka}`, basal sliding ramped up during the first 5 ka

       * modified, low-resolution version of the spin-up for ISMIP6 InitMIP (Greve et al. :cite:`greve_etal_2017a`)

     * ``v5_ant64_b2_spinup09_init100a``, 
       ``v5_ant64_b2_spinup09_fixtopo``, ``v5_ant64_b2_spinup09`` and ``v5_ant64_b2_future09_ctrl``

       * Antarctic ice sheet with hybrid shallow-ice--shelfy-stream dynamics (Bernales et al. :cite:`bernales_etal_2017a`) and ice shelves (SSA)

       * Resolution 64 km, :math:`t=-140.1\ldots{}-140.0\,\mathrm{ka}` for the init run without basal sliding (..._init100a)

       * :math:`t=-140\ldots{}0\,\mathrm{ka}` for the run with almost fixed topography (..._fixtopo)

       * basal sliding ramped up during the first 5 ka

       * :math:`t=-0.5\ldots{}0\,\mathrm{ka}` for the final, freely-evolving-topography part of the (..._spinup09),

       * :math:`t=0\ldots{}100\,\mathrm{a}` for the constant-climate control run (..._future09_ctrl)

       * 64-km version of the spin-up and the constant-climate control run for ISMIP6 InitMIP; Greve and Galton-Fenzi (pers. comm. 2017).

     * ``v5_asf2_steady and v5_asf2_surge``

       * Austfonna, SIA, resolution 2 km, :math:`t=0\ldots{}10\,\mathrm{ka}`

       * Similar to Dunse et al. :cite:`dunse_etal_2011`'s Exp. 2 (steady fast flow) and Exp. 5 (surging-type flow), respectively      

     * ``v5_nmars10_steady``, ``v5_smars10_steady``

       * North-/south-polar cap of Mars, SIA, resolution 10 km, :math:`t=-10\,\mathrm{Ma}\ldots{}0`

       * Steady-state runs by Greve :cite:`greve_2007b`
 
     * ``v5_nhem80_nt012_new``

       * northern hemisphere, SIA, resolution 80 km, :math:`t=-250\ldots{}0\,\mathrm{ka}`

       * Similar to run nt012 by Greve et al. :cite:`greve_etal_1999a`

     * ``v5_heino50_st``

       * ISMIP HEINO standard run ST

       * SIA, resolution 50 km, :math:`t=0\ldots{}200\,\mathrm{ka}` (Calov et al. :cite:`calov_etal_2010`).
