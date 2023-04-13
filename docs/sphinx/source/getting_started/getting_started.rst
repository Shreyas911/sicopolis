.. _getting_started:

Getting started
***************

Requirements
============

Lorem ipsum...

Download
========

Lorem ipsum...

Initial configuration
=====================

Lorem ipsum...

Directory structure
===================

* ``runs`` :

  1. Configuration file ``sico_configs.sh``.

  2. Shell script (bash) ``sico.sh`` for running a single simulation.

  3. Shell scripts (bash) ``multi_sico_1.sh`` and ``multi_sico_2.sh`` for running multiple simulations by repeated calls of ``sico.sh``.

  4. Subdirectory ``headers``: specification files ``sico_specs_{run_name}.h`` ({run_name}: name of run) for a number of computationally rather inexpensive test runs.

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

       * Resolution 64 km, :math:`t=-140.1\ldots{}-140.0\,\mathrm{ka}` for the init run without basal sliding (...\_init100a)

       * :math:`t=-140\ldots{}0\,\mathrm{ka}` for the run with almost fixed topography (...\_fixtopo)

       * basal sliding ramped up during the first 5 ka

       * :math:`t=-0.5\ldots{}0\,\mathrm{ka}` for the final, freely-evolving-topography part of the (...\_spinup09),

       * :math:`t=0\ldots{}100\,\mathrm{a}` for the constant-climate control run (...\_future09\_ctrl)

       * 64-km version of the spin-up and the constant-climate control run for ISMIP6 InitMIP; Greve and Galton-Fenzi (pers.\ comm.\ 2017).

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

* ``src`` :

  * Directory that contains the main program file sicopolis.F90.

  1. Subdirectory ``subroutines/general`` : general subroutines, for any modelled domain.
  
  2. Subdirectory ``subroutines/ant`` : subroutines specific for the Antarctic ice sheet.

  3. Subdirectory ``subroutines/emtp2sge`` : subroutines specific for the EISMINT Phase 2 Simplified Geometry Experiments.

  4. Subdirectory ``subroutines/grl`` : subroutines specific for the Greenland ice sheet.

  5. Accordingly subdirectories subroutines/asf, nhem, scand, tibet, nmars and smars for Austfonna, the northern hemisphere, Scandinavia, Tibet and the north and south polar caps of Mars, respectively.

  6. Subdirectory ``subroutines/tapenade`` : AD specific subroutines and files.

  7. Subdirectory ``subroutines/xyz`` : Framework to create new domains, this directory is empty by default.

* ``sico_in`` :

  * Directory that contains input data files for SICOPOLIS.

  1. Subdirectory ``general`` : general input files, for any modelled domain.

  2. Subdirectory ``ant`` : input files specific for the Antarctic ice sheet. 

  3. Subdirectory ``emtp2sge`` : input files specific for the EISMINT Phase 2 Simplified Geometry Experiments.

  4. Subdirectory ``grl``: input files specific for the Greenland ice sheet.

  5. Accordingly subdirectories asf, nhem, scand, tibet, nmars and smars for Austfonna, the northern hemisphere, Scandinavia, Tibet and the north and south polar caps of Mars, respectively.

  6. Subdirectory \textbf{xyz}: Framework to create new domains, place your input files here.

* ``test_ad`` :

  * AD specific utilities and CI testing framework

* ``sico_out`` :

  * Empty directory into which output files of SICOPOLIS simulations are written.

* ``docs`` :

  * Documentation with quick-start manual, sphinx docs, JOSS paper, doxygen, etc.

* ``tools`` :

  * Tools to help with forward modeling, eg - ``resolution_doubler`` , ``make_ismip_output`` , etc.

  1. Program ``make_ismip_output``

     * Generating ISMIP output (see http://tinyurl.com/clic-ismip6) from the NetCDF time-slice files produced by SICOPOLIS

     * For simulation run ``./tools.sh -p make_ismip_output -m run_name``

     * For further options, try ``./tools.sh -h``

  2. Program ``resolution_doubler``

     * Doubling the horizontal resolution of a NetCDF time-slice output file produced by SICOPOLIS

     * For simulation run name, to be executed by ``./tools.sh -p resolution_doubler -m run_name`` 

     * For further options, try ``./tools.sh -h`` 

     * For example, run ``v5_grl10_b2_paleo21`` (10 km resolution) requires the resolution doubled output of run ``v5_grl20_b2_paleo21`` (20 km resolution) for :math:`t=-9 \mathrm{ka}` as initial condition. In order to create it, execute the resolution doubler for run ``v5_grl20_b2_paleo21`` (i.e., with the option ``-m v5_grl20_b2_paleo21``) and enter 

       * Number of time-slice file (with leading zeros, 4 digits) :math:`> 0004` 

     * This will convert the original time-slice file ``v5_grl20_b2_paleo210004.nc`` to the resolution-doubled file ``v5_grl20_b2_paleo21_dbl_0004.nc`` that serves as initial conditions for run ``v5_grl10_b2_paleo21`` .

How to run a simulation
=======================

Lorem ipsum...

Output files
============

Lorem ipsum...

Plotting
========

Lorem ipsum...

Some useful tools
=================

Lorem ipsum...
