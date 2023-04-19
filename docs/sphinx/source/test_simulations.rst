.. _test_simulations:

Test simulations
****************

These are a number of computationally rather inexpensive test simulations, of which the :ref:`run-specs header files <getting_started-run_specs_headers>` are contained in the SICOPOLIS repository.

Run v5_vialov3d25
  | 3-d version of the 2-d Vialov profile (Vialov :cite:`vialov_1958`),
  | SIA, resolution 25 km, :math:`t=0\ldots{}100\,\mathrm{ka}`.
  | Similar to the EISMINT Phase 1 fixed-margin experiment (Huybrechts et al. :cite:`huybrechts_etal_1996`), but without thermodynamics. Instead, isothermal conditions with :math:`T=-10^{\circ}\mathrm{C}` everywhere are assumed.

Run v5_emtp2sge25_expA
  | EISMINT Phase 2 Simplified Geometry Experiment A,
  | SIA, resolution 25 km, :math:`t=0\ldots{}200\,\mathrm{ka}` (Payne et al. :cite:`payne_etal_2000`).
  | The thermodynamics solver for this run is the one-layer melting-CTS enthalpy scheme (ENTM), while all other runs employ the polythermal two-layer scheme (POLY) (Greve and Blatter :cite:`greve_blatter_2016`).

Run v5_grl16_bm5_ss25ka
  | Greenland ice sheet, SIA, resolution 16 km,
  | short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for modern climate conditions (unpublished).

Run v5_ant40_b2_ss25ka
  | Antarctic ice sheet without ice shelves, SIA, resolution 40 km,
  | short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for modern climate conditions (unpublished).

Run v5_grl20_b2_paleo21
  | Greenland ice sheet, SIA, resolution 20 km,
  | :math:`t=-140\,\mathrm{ka}\ldots{}0`, basal sliding ramped up during the first 5 ka.
  | Modified, low-resolution version of the spin-up for ISMIP6 InitMIP (Greve et al. :cite:`greve_etal_2017a`).

Runs v5\_grl10\_b2\_\{paleo21, future21\_ctrl, future21\_asmb\}
  | Greenland ice sheet, SIA, resolution 10 km,
  | :math:`t=-9\,\mathrm{ka}\ldots{}0` for the paleo run, :math:`t=0\ldots{}100\,\mathrm{a}` for the future runs.
  | 10-km version of the spin-up and the schematic future climate runs for ISMIP6 InitMIP (Greve et al. :cite:`greve_etal_2017a`).

Runs v5\_ant64\_b2\_\{spinup09\_init100a, spinup09\_fixtopo, spinup09, future09\_ctrl\}
  | Antarctic ice sheet with hybrid shallow-ice--shelfy-stream dynamics
  | (Bernales et al. :cite:`bernales_etal_2017a`) and ice shelves (SSA), resolution 64 km;
  | :math:`t=-140.1\ldots{}-140.0\,\mathrm{ka}` for the init run without basal sliding (..._init100a),
  | :math:`t=-140\,\mathrm{ka}\ldots{}0` for the run with almost fixed topography (..._fixtopo), basal sliding ramped up during the first 5 ka,
  | :math:`t=-0.5\,\mathrm{ka}\ldots{}0` for the final, freely-evolving-topography part of the (..._spinup09),
  | :math:`t=0\ldots{}100\,\mathrm{a}` for the constant-climate control run (..._future09_ctrl).
  | 64-km version of the spin-up and the constant-climate control run for ISMIP6 InitMIP (Greve and Galton-Fenzi, pers. comm. 2017).

Runs v5_asf2_steady, v5_asf2_surge
  | Austfonna, SIA, resolution 2 km, :math:`t=0\ldots{}10\,\mathrm{ka}`.
  | Similar to Dunse et al. :cite:`dunse_etal_2011`'s Exp. 2 (steady fast flow) and Exp. 5 (surging-type flow), respectively.

Runs v5_nmars10_steady, v5_smars10_steady
  | North-/south-polar cap of Mars, SIA, resolution 10 km, :math:`t=-10\,\mathrm{Ma}\ldots{}0`.
  | Steady-state runs by Greve :cite:`greve_2007b`.
 
Run v5_nhem80_nt012_new
  | Northern hemisphere, SIA, resolution 80 km, :math:`t=-250\,\mathrm{ka}\ldots{}0`.
  | Similar to run nt012 by Greve et al. :cite:`greve_etal_1999a`.

Run v5_heino50_st
  | ISMIP HEINO standard run ST, SIA, resolution 50 km, :math:`t=0\ldots{}200\,\mathrm{ka}` (Calov et al. :cite:`calov_etal_2010`).

-------------

**Computing times:**

+-----------------------------------+------------+---------------------+--------------------+
| Run                               | Model time | Time step\ :sup:`†` | CPU time\ :sup:`‡` |
+===================================+============+=====================+====================+
| v5\_vialov3d25                    | 100 ka     | 20 a                | 1.0 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_emtp2sge25\_expA              | 200 ka     | 20 a                | 3.9 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_grl16\_bm5\_ss25ka            | 25 ka      | 5 a                 | 9.7 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_ant40\_b2\_ss25ka             | 25 ka      | 10 a                | 5.0 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_grl20\_b2\_paleo21            | 140 ka     | 5 a                 | 0.8 hrs            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_grl10\_b2\_paleo21\ :sup:`\*` | 9 ka       | 1 a                 | 1.0 hrs            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_grl10\_b2\_future21\_ctrl     | 100 a      | 1 a                 | 0.9 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_grl10\_b2\_future21\_asmb     | 100 a      | 1 a                 | 0.9 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_ant64\_b2\_spinup09\_init100a | 100 a      | 2 / 10 a\ :sup:`†`  | 4.1 sec            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_ant64\_b2\_spinup09\_fixtopo  | 140 ka     | 5 / 10 a\ :sup:`†`  | 0.7 hrs            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_ant64\_b2\_spinup09           | 500 a      | 2 / 10 a\ :sup:`†`  | 0.5 min            |
+-----------------------------------+------------+---------------------+--------------------+
| v5\_ant64\_b2\_future09\_ctrl     | 100 a      | 2 / 10 a\ :sup:`†`  | 6.1 sec            |
+-----------------------------------+------------+---------------------+--------------------+

| Table 1: Model times, time steps and computing (CPU) times for the EISMINT, Greenland and Antarctica test simulations contained in the script multi_sico_1.sh, run with SICOPOLIS V5-dev (branch develop, revision 9c909c3c2) and the Intel Fortran Compiler 19.1 for Linux (optimization options -xHOST -O3 -no-prec-div) on a 12-Core Intel Xeon Gold 6256 (3.6 GHz) PC under openSUSE Leap 15.4.
| \ :sup:`†`: If one value is given, this is the common dynamic (velocity, ice thickness) and thermodynamic (temperature, water content, age) time step. If two values are given (marked by the dagger (\ :sup:`†`) symbol), the first one is the dynamic, the second one the thermodynamic time step.
| \ :sup:`‡`: All runs were done on a single core only. The v5_ant64_b2_xxx runs that include ice shelves can be done on multiple cores using OpenMP for the SSA solver. However, at the employed, low resolution of 64 km the solver does not scale well, and the gain in wall clock time by using multiple cores is very small.
| \ :sup:`\*`: For this run, see the remark in the :ref:`subsection on the resolution-doubler tool <plotting_and_tools-res_dbl>`.
