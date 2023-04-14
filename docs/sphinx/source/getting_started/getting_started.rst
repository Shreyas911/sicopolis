.. _getting_started:

Getting started
***************

Requirements
============

Unix-like system (e.g., Linux).

| Fortran compiler.
| So far, the GNU GCC (gfortran) and Intel Fortran (ifort) compilers are supported. If you wish to use a different compiler, please contact help@sicopolis.net.

SICOPOLIS writes output in `NetCDF <https://doi.org/10.5065/D6H70CW6>`__ format (plus some ASCII). An installation of NetCDF version 3.6.x or newer is therefore required. For installation tips, see :ref:`Dependencies/NetCDF <dependencies-netcdf>`.

For the shallow-shelf/shelfy-stream solver, the Library of Iterative Solvers for Linear Systems (`Lis <https://www.ssisc.org/lis/>`__, version 1.4.43 or newer) is required. For installation tips, see :ref:`Dependencies/Lis <dependencies-lis>`.

Download
========

Clone the Git repository (recommended)
  The Git repository of SICOPOLIS is kindly hosted by the GitLab system of the Alfred Wegener Institute for Polar and Marine Research (AWI) in Bremerhaven, Germany. Front page: https://gitlab.awi.de/sicopolis/sicopolis/.

  Cloning the latest develop revision::

    git clone --branch develop
        https://gitlab.awi.de/sicopolis/sicopolis.git

  (Cloning with SSH instead of HTTPS is also available. See the above GitLab front page link for details.)

  You should then have a new directory "sicopolis" that contains the entire program package.

Download an archive
  Go to https://gitlab.awi.de/sicopolis/sicopolis/-/tree/develop, click on the download symbol and choose the desired format (zip, tar.gz, tar.bz2, tar).

  Tagged versions of SICOPOLIS are also available from `Zenodo <https://doi.org/10.5281/zenodo.3687337>`__.

Initial configuration
=====================

1. Go to the new directory "sicopolis" and execute the following bash scripts::

      ./copy_templates.sh
      ./get_input_files.sh

   The latter can be configured if you want to download only selected input files (default is downloading everything). To do so, open it with a text editor and change the flag variables before execution.

2. Locate the file sico_configs.sh in the directory sicopolis/runs, and open it with a text editor.

   Set the flags LIS_FLAG, OPENMP_FLAG and LARGE_DATA_FLAG according to your needs. 

   Default is "true"/"true"/"false", which works for all test simulations included in the SICOPOLIS package. LIS_FLAG and OPENMP_FLAG can be set to "false" for simulations with pure shallow-ice dynamics. However, "true" is required for simulations with shallow-shelf dynamics (for floating ice) or hybrid shallow-ice--shelfy-stream dynamics (for grounded ice). For high-resolution simulations (e.g., Greenland/5 km or Antarctica/8 km), LARGE_DATA_FLAG must be set to "true".

   Set NETCDFHOME to the correct path of your NetCDF installation.

   If LIS_FLAG = "true", set LISHOME to the correct path of your Lis installation.

   Depending on your system, some additional settings might have to be added in sico_configs.sh ("module load" commands for dynamic loading etc.).

3. Locate the file sico_environment.sh in the directory sicopolis/runs, open it with a text editor, and replace the "Default" entry for SICO_INSTITUTION by the name of your institution (max. 256 characters).

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

  * Documentation for SICOPOLIS.

* ``tools`` :

  * Directory that contains some useful tools and a shell script (tools.sh) to execute them (see ":ref:`plotting_and_tools`").

How to run a simulation
=======================

Lorem ipsum...

.. _getting_started-output:

Output files
============

Lorem ipsum...
