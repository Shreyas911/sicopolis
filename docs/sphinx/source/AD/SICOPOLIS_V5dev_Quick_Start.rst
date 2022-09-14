--------------

| **SICOPOLIS V5-dev
  – Quick Start Manual –**

--------------

Ralf Greve

--------------

| Institute of Low Temperature Science, Hokkaido University,
| Kita-19, Nishi-8, Kita-ku, Sapporo 060-0819, Japan

--------------

--------------

| Copyright 2009–2022 Ralf Greve
| (with contributions by Jorge Bernales, Sebastian Beyer, Heinz Blatter,
  Reinhard Calov, Thorben Dunse, Eduardo Flandez, Ben Galton-Fenzi,
  Thomas Gölles, Björn Grieger, Philipp Hancke, Patrick Heimbach, Nina
  Kirchner, Thomas Kleiner, Sascha Knell, Anne Le Brocq, Liz Curry
  Logan, Sri Hari Krishna Narayanan, Alex Robinson, Fuyuki Saito,
  Tatsuru Sato, Marius Schäfer, Matthias Scheiter, Oliver J. Stenzel,
  Malte Thoma, Roland Warner)

This file is part of SICOPOLIS.

SICOPOLIS is free software. It can be redistributed and/or modified
under the terms of the GNU General Public License
(http://www.gnu.org/licenses/) as published by the Free Software
Foundation, either version 3 of the License, or (at the user’s option)
any later version.

SICOPOLIS is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

--------------

Requirements
============

-  Unix-like system (e.g., Linux).

-  | Fortran compiler.
   | So far, the GNU GCC (gfortran) and Intel Fortran (ifort) compilers
     are supported. If you wish to use a different compiler, please
     contact :math:`<`\ help@sicopolis.net\ :math:`>`.

-  SICOPOLIS writes output in NetCDF format (plus some ASCII). An
   installation of NetCDF version 3.6.x or newer
   (https://www.unidata.ucar.edu/software/netcdf/) is therefore
   required. For installation support, see
   Appendix `9 <#sect_install_nc>`__.

-  For the shallow-shelf/shelfy-stream solver, the Library of Iterative
   Solvers for Linear Systems (Lis, version 1.4.43 or newer) is required
   (https://www.ssisc.org/lis/). For installation support, see
   Appendix `10 <#sect_install_lis>`__.

Download
========

-  **Using Git**

   The Git repository of SICOPOLIS is kindly hosted by the GitLab system
   of the Alfred Wegener Institute for Polar and Marine Research (AWI)
   in Bremerhaven, Germany. Front page:
   https://gitlab.awi.de/sicopolis/sicopolis/ .

   Cloning the latest develop revision:

   | ``git clone --branch develop \``
   | ``https://gitlab.awi.de/sicopolis/sicopolis.git``

   (Cloning with SSH instead of HTTPS is also available. See the above
   GitLab front page link for details.)

   You should then have a new directory “sicopolis” that contains the
   entire program package.

Initial configuration
=====================

#. Go to the new directory “sicopolis” and execute the following bash
   scripts:

   | ``./copy_templates.sh``
   | ``./get_input_files.sh``

   The latter can be configured if you want to download only selected
   input files (default is downloading everything). To do so, open it
   with a text editor and change the flag variables before execution.

#. Locate the file sico_configs.sh in the directory sicopolis/runs, and
   open it with a text editor.

#. | Set the flags
   | LIS_FLAG, OPENMP_FLAG and LARGE_DATA_FLAG
   | according to your needs.

   Default is “true”/“true”/“false”, which works for all test
   simulations included in the SICOPOLIS package. LIS_FLAG and
   OPENMP_FLAG can be set to “false” for simulations with pure
   shallow-ice dynamics. However, “true” is required for simulations
   with shallow-shelf dynamics (for floating ice) or hybrid
   shallow-ice–shelfy-stream dynamics (for grounded ice). For
   high-resolution simulations (e.g., Greenland/5 km or
   Antarctica/8 km), LARGE_DATA_FLAG must be set to “true”.

#. Set NETCDFHOME to the correct path of your NetCDF installation.

   If LIS_FLAG = “true”, set LISHOME to the correct path of your Lis
   installation.

#. Depending on your system, some additional settings might have to be
   added in sico_configs.sh (``module load``  commands for dynamic
   loading etc.).

#. Locate the file sico_environment.sh in the directory sicopolis/runs,
   open it with a text editor, and replace the “Default” entry for
   SICO_INSTITUTION by the name of your institution (max. 256
   characters).

Files and directories in “sicopolis”
====================================

-  **runs**:

   Configuration file sico_configs.sh.

   Shell script (bash) sico.sh for running a single simulation.

   Shell scripts (bash) multi_sico_1.sh and multi_sico_2.sh for running
   multiple simulations by repeated calls of sico.sh.

   Subdirectory **headers**: specification files
   sico_specs\_\ *run_name*.h (*run_name*: name of run) for a number of
   computationally rather inexpensive test runs.

   -  | Run v5_vialov3d25
      | :math:`\longrightarrow` 3-d version of the 2-d “Vialov profile”
        :raw-latex:`\citep{vialov_1958}`,
      | SIA, resolution 25 km, :math:`t=0\ldots{}100\,\mathrm{ka}`.
      | Similar to the EISMINT Phase 1 fixed-margin experiment
      | :raw-latex:`\citep{huybrechts_etal_1996}`, but without
        thermodynamics. Instead,
      | isothermal conditions with
        :math:`T=-10\ensuremath{^\circ\mathrm{C}}` everywhere are
        assumed.

   -  | Run v5_emtp2sge25_expA
      | :math:`\longrightarrow` EISMINT Phase 2 Simplified Geometry
        Experiment A,
      | SIA, resolution 25 km, :math:`t=0\ldots{}200\,\mathrm{ka}`
        :raw-latex:`\citep{payne_etal_2000}`.
      | The thermodynamics solver for this run is the one-layer
        melting-CTS
      | enthalpy scheme (ENTM), while all other runs employ the
      | polythermal two-layer scheme (POLY)
        :raw-latex:`\citep{greve_blatter_2016}`.

   -  | Run v5_grl16_bm4_ss25ka
      | :math:`\longrightarrow` Greenland ice sheet, SIA, resolution
        16 km,
      | short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for
        modern climate conditions
      | (unpublished).

   -  | Run v5_ant40_b2_ss25ka
      | :math:`\longrightarrow` Antarctic ice sheet without ice shelves,
        SIA, resolution 40 km,
      | short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for
        modern climate conditions
      | (unpublished).

   -  | Run v5_grl20_b2_paleo21
      | :math:`\longrightarrow` Greenland ice sheet, SIA, resolution
        20 km,
      | :math:`t=-140\ldots{}0\,\mathrm{ka}`, basal sliding ramped up
        during the first 5 ka
      | [modified, low-resolution version of the spin-up for ISMIP6
        InitMIP;
      | :raw-latex:`\citet{greve_etal_2017a}`].

   -  | Runs v5_grl10_b2_paleo21 and
      | v5_grl10_b2_future21_ctrl/..._asmb
      | :math:`\longrightarrow` Greenland ice sheet, SIA, resolution
        10 km,
      | :math:`t=-9\ldots{}0\,\mathrm{ka}` for the paleo run,
        :math:`t=0\ldots{}100\,\mathrm{a}` for the future runs
      | [10-km version of the spin-up and the schematic future climate
        runs for
      | ISMIP6 InitMIP; :raw-latex:`\citet{greve_etal_2017a}`].

   -  | Runs v5_ant64_b2_spinup09_init100a,
      | v5_ant64_b2_spinup09_fixtopo, v5_ant64_b2_spinup09 and
      | v5_ant64_b2_future09_ctrl
      | :math:`\longrightarrow` Antarctic ice sheet with hybrid
        shallow-ice–shelfy-stream dynamics
      | :raw-latex:`\citep{bernales_etal_2017}` and ice shelves (SSA),
        resolution 64 km,
      | :math:`t=-140.1\ldots{}-140.0\,\mathrm{ka}` for the init run
        without basal sliding (..._init100a),
      | :math:`t=-140\ldots{}0\,\mathrm{ka}` for the run with almost
        fixed topography (..._fixtopo),
      | basal sliding ramped up during the first 5 ka,
      | :math:`t=-0.5\ldots{}0\,\mathrm{ka}` for the final,
        freely-evolving-topography part of the
      | spin-up (..._spinup09),
      | :math:`t=0\ldots{}100\,\mathrm{a}` for the constant-climate
        control run (..._future09_ctrl)
      | [64-km version of the spin-up and the constant-climate control
        run for
      | ISMIP6 InitMIP; Greve and Galton-Fenzi (pers. comm. 2017)].

   -  | Runs v5_asf2_steady and v5_asf2_surge
      | :math:`\longrightarrow` Austfonna, SIA, resolution 2 km,
        :math:`t=0\ldots{}10\,\mathrm{ka}`
      | [similar to :raw-latex:`\citeauthor{dunse_etal_2011}`’s
        (:raw-latex:`\citeyear{dunse_etal_2011}`) Exp. 2 (steady fast
        flow) and
      | Exp. 5 (surging-type flow), respectively].

   -  | Runs v5_nmars10_steady and v5_smars10_steady
      | :math:`\longrightarrow` North-/south-polar cap of Mars, SIA,
        resolution 10 km, :math:`t=-10\,\mathrm{Ma}\ldots{}0`
      | [steady-state runs by :raw-latex:`\citet{greve_2007b}`].

   -  | Run v5_nhem80_nt012_new
      | :math:`\longrightarrow` northern hemisphere, SIA, resolution
        80 km, :math:`t=-250\ldots{}0\,\mathrm{ka}`
      | [similar to run nt012 by :raw-latex:`\citet{greve_etal_1999a}`].

   -  | Run v5_heino50_st
      | :math:`\longrightarrow` ISMIP HEINO standard run ST,
      | SIA, resolution 50 km, :math:`t=0\ldots{}200\,\mathrm{ka}`
        :raw-latex:`\citep{calov_etal_2010}`.

-  **src**:

   Directory that contains the main program file sicopolis.F90.

   -  Subdirectory **subroutines/general**: general subroutines, for any
      modelled domain.

   -  Subdirectory **subroutines/ant**: subroutines specific for the
      Antarctic ice sheet.

   -  Subdirectory **subroutines/emtp2sge**: subroutines specific for
      the EISMINT Phase 2 Simplified Geometry Experiments.

   -  Subdirectory **subroutines/grl**: subroutines specific for the
      Greenland ice sheet.

   -  Accordingly subdirectories subroutines/asf, nhem, scand, tibet,
      nmars and smars for Austfonna, the northern hemisphere,
      Scandinavia, Tibet and the north and south polar caps of Mars,
      respectively.

   -  Subdirectory **subroutines/xyz**: see Appendix `11 <#sect_xyz>`__.

-  **sico_in**:

   Directory that contains input data files for SICOPOLIS.

   -  Subdirectory **general**: general input files, for any modelled
      domain.

   -  Subdirectory **ant**: input files specific for the Antarctic ice
      sheet.

   -  Subdirectory **emtp2sge**: input files specific for the EISMINT
      Phase 2 Simplified Geometry Experiments.

   -  Subdirectory **grl**: input files specific for the Greenland ice
      sheet.

   -  Accordingly subdirectories asf, nhem, scand, tibet, nmars and
      smars for Austfonna, the northern hemisphere, Scandinavia, Tibet
      and the north and south polar caps of Mars, respectively.

   -  Subdirectory **xyz**: see Appendix `11 <#sect_xyz>`__.

-  **sico_out**:

   Empty directory into which output files of SICOPOLIS simulations are
   written.

-  **docu**:

   Directory that contains some documentation.

   -  | Subdirectory **quick_start**:
      | LaTeX source for this manual (PDF must be built with make).

   -  | Subdirectory **doxygen**: documentation to be created by Doxygen
      | (optional, see doxygen-config/README.md).

      -  html/index.html :math:`\longrightarrow` Source code browser.

      -  latex/refman.pdf :math:`\longrightarrow` Reference manual.

-  **tools**:

   See Sect. `8 <#sect_tools>`__.

-  **license**:

   Directory that contains a copy of the GNU General Public License
   (version 3).

How to run a simulation
=======================

#. | For example, to run simulation v5_grl16_bm4_ss25ka, use the script
     sico.sh:
   | ``(./sico.sh -m v5_grl16_bm4_ss25ka) >out_001.dat 2>&1 &``
   | (from directory sicopolis/runs, bash required). Accordingly for the
     other simulations.
   | For further options, try ``./sico.sh -h``.
   | :math:`\backslash`!/ Do not use out\_\ *run_name*.dat for the
     redirected output of sico.sh.
   | This name is reserved for the runtime output of SICOPOLIS itself.
   | (Both are very useful in case of compilation or runtime errors!)

#. | Alternatively, if you prefer to run all simulations consecutively,
     execute the script multi_sico_1.sh:
   | ``(./multi_sico_1.sh) >out_multi_100.dat 2>&1 &``
   | For further options, try ``./multi_sico_1.sh -h``.

Computing times
~~~~~~~~~~~~~~~

The approximate computing times for the simulations are listed in
Table `1 <#table_comp_times>`__ (Appendix `12 <#sect_table>`__).

.. _sect_output:

Output files
============

Output files are written by default to the directory
sicopolis/sico_out/*run_name* (this can be changed with the ``-d``
option). Four types are produced:

-  **run_name.log**:

   ASCII file that lists the main specifications of simulation
   *run_name*.

-  **run_name.ser, run_name_ser.nc**:

   Time-series files (ASCII, NetCDF) that contain scalar variables:

   -  Time, t

   -  Surface temperature anomaly, D_Ts, or glacial index, glac_ind
      (forcing)

   -  Sea level, z_sl (forcing)

   -  Total ice volume, V

   -  Volume of grounded ice, V_g

   -  Volume of floating ice, V_f

   -  Total ice area, A

   -  Area of grounded ice, A_g

   -  Area of floating ice, A_f

   -  Ice volume above flotation in sea level equivalent, V_sle

   -  Volume of temperate ice, V_t

   -  Area of temperate-based grounded ice, A_t

   -  Maximum ice thickness, H_max

   -  Maximum thickness of temperate ice, H_t_max

   -  Maximum surface elevation, zs_max

   -  Maximum surface speed, vs_max

   -  Maximum basal temperature (relative to pmp), Tbh_max

   -  (Some more in the NetCDF file, try ``ncdump -h run_name_ser.nc``)

-  **run_name.core, run_name_core.nc**:

   Time-series files (ASCII, NetCDF) that contain for selected locations
   xxx:

   -  Time, t

   -  Surface temperature anomaly, D_Ts, or glacial index, glac_ind
      (forcing)

   -  Sea level, z_sl (forcing)

   -  Thickness, H_xxx

   -  Surface velocity, v_xxx

   -  Basal temperature, T_xxx

   -  (Some more in the NetCDF file, try ``ncdump -h run_name_core.nc``)

   | For the Greenland ice sheet, these data are written for seven
     locations:
   | GRIP (xxx=GR), GISP2 (xxx=G2), Dye 3 (xxx=D3), Camp Century
     (xxx=CC), NorthGRIP (xxx=NG), NEEM (xxx=NE), EastGRIP (xxx=EG).

   | For the Antarctic ice sheet, these data are written for six
     locations:
   | Vostok (xxx=Vo), Dome A (xxx=DA), Dome C (xxx=DC), Dome F (xxx=DF),
     Kohnen (xxx=Ko), Byrd (xxx=By).

-  **run_name0001.nc**, **run_name0002.nc**, ...:

   Complete set of fields (topography, velocity, temperature etc.,
   written in NetCDF (*.nc) format) for selected time slices defined in
   specifications file.

   For example, simulation v5_grl16_bm4_ss25ka produces three files
   v5_grl16_bm4_ss25ka0001.nc, v5_grl16_bm4_ss25ka0002.nc and
   v5_grl16_bm4_ss25ka0003.nc, which correspond to :math:`t=0`, 10 ka
   and 25 ka, respectively.

.. _sect_plotting:

Plotting
========

The output described in Sect. `6 <#sect_output>`__ can be visualized
with any plotting tool at the user’s preference. Ncview
(http://meteora.ucsd.edu/:math:`\sim`\ pierce/ncview_home_page.html) is
a very nice browser for NetCDF files to get a quick and easy look. For
more sophisticated plots, one possibility is to use MATLAB, which has an
extensive library for NetCDF files
(https://www.mathworks.com/help/matlab/network-common-data-form.html).
For instance, the following script plots the final surface topography of
the Greenland simulation v5_grl16_bm4_ss25ka (credit: Mathieu Morlighem,
University of California Irvine).

| ``filename = 'v5_grl16_bm4_ss25ka0003.nc';``
| ``x = ncread(filename,'x');``
| ``y = ncread(filename,'y');``
| ``surf = ncread(filename,'zs');``
| ``% Display surface elevation``
| ``%    (transposition needed because MATLAB is column-oriented)``
| ``imagesc(x*1e-3,y*1e-3,surf'); axis xy equal; caxis([0 3200]); colorbar``

.. _sect_tools:

Some useful tools
=================

The directory sicopolis/tools contains some useful tools.

.. _ssect_make_ismip_output:

Program make_ismip_output
-------------------------

| Generating ISMIP output (see http://tinyurl.com/clic-ismip6) from the
  NetCDF time-slice files produced by SICOPOLIS (see
  Sect. `6 <#sect_output>`__). For simulation *run_name*, to be executed
  by
| ``./tools.sh -p make_ismip_output -m run_name``
| For further options, try ``./tools.sh -h``.

.. _ssect_resolution_doubler:

Program resolution_doubler
--------------------------

| Doubling the horizontal resolution of a NetCDF time-slice output file
  produced by SICOPOLIS (see Sect. `6 <#sect_output>`__). For simulation
  *run_name*, to be executed by
| ``./tools.sh -p resolution_doubler -m run_name``
| For further options, try ``./tools.sh -h``.

| For example, run v5_grl10_b2_paleo21 (10 km resolution) requires the
  resolution-doubled output of run v5_grl20_b2_paleo21 (20 km
  resolution) for :math:`t=-9\,\mathrm{ka}` as initial condition. In
  order to create it, execute the resolution doubler for run
  v5_grl20_b2_paleo21 (i.e., with the option ``-m v5_grl20_b2_paleo21``)
  and enter
| Number of time-slice file (with leading zeros, 4 digits) :math:`>`
  ``0004``
| This will convert the original time-slice file
  v5_grl20_b2_paleo210004.nc to the resolution-doubled file
  v5_grl20_b2_paleo21_dbl_0004.nc that serves as initial conditions for
  run v5_grl10_b2_paleo21.

.. _sect_install_nc:

Installation of NetCDF
======================

NetCDF (Network Common Data Form) is a common format for scientific data
(https://www.unidata.ucar.edu/software/netcdf/) that is also used by
SICOPOLIS. The NetCDF C and Fortran libraries are required.

For **GCC**, installation from a package manager is recommended. Under
openSUSE Leap 15.2, install netcdf, netcdf-devel, netcdf-devel-static,
netcdf-fortran, netcdf-fortran-devel, netcdf-fortran-static, ncview.
This requires the repositories “Software for Scientists and Engineers”
and “sebschub’s Home Project”. Details (especially the required
repositories) will differ for other systems.

For the **Intel compiler**, manual installation is required. The C and
Fortran libraries are available for download on the NetCDF website as
zip or tar archives. Unzip them into temporary source directories.

-  Prior to version 4.2, a single archive contained both the C and
   Fortran libraries. A minimal installation for version 4.1.3 (without
   NetCDF-4 support) can be done by changing to the source directory,
   then:

   ::

        export NCDIR=/opt/netcdf
        export CC=icc
        export FC=ifort
        export CFLAGS="-O2"
        export CPPFLAGS=
        export FCFLAGS="-O2"
        export FFLAGS=${FCFLAGS}
        ./configure --prefix=${NCDIR} --libdir=${NCDIR}/lib \
                    --disable-netcdf-4
        make install

-  Since version 4.2, the C and Fortran libraries must be installed
   separately. If NetCDF-4 support is dispensable, the following
   installation should work (tested under openSUSE Leap 15.2 and
   icc/ifort 19.1 with versions netcdf-c-4.8.0 and netcdf-fortran-4.5.3
   as of January 25, 2021).

   Step 1: Change to the source directory of the C library, then:

   ::

        export NCDIR=/opt/netcdf
        export CC=icc
        export FC=ifort
        export CFLAGS="-O2"
        export CPPFLAGS=
        export FCFLAGS="-O2"
        export FFLAGS=${FCFLAGS}
        ./configure --prefix=${NCDIR} --libdir=${NCDIR}/lib \
                    --disable-netcdf-4 --enable-logging
        make install

   Step 2: Change to the source directory of the Fortran library, then:

   ::

        export NFDIR=/opt/netcdf
        export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
        export CPPFLAGS=-I${NCDIR}/include
        export LDFLAGS=-L${NCDIR}/lib
        ./configure --prefix=${NFDIR} --libdir=${NFDIR}/lib \
                    --disable-netcdf-4 --enable-logging
        make install

-  For a complete build with NetCDF-4 support, additional libraries are
   required. See the NetCDF website for further instructions.

Installation under /opt usually requires admin rights. The same holds
for the common alternative /usr/local. For a local installation, replace
it by ‘/home/:math:`<`\ my_user_name\ :math:`>`/local’.

.. _sect_install_lis:

Installation of Lis
===================

Lis (Library of Iterative Solvers for linear systems) is a software
library for solving discretized linear equations
:raw-latex:`\citep{nishida_2010}`.

Download the latest version of Lis as a zip archive from
https://www.ssisc.org/lis/ (as of January 23, 2021: lis-2.0.30.zip).
Unzip the archive into a temporary directory.

For **GCC**, install lis by executing:

::

     export LISDIR=/opt/lis
     ./configure --prefix=${LISDIR} --libdir=${LISDIR}/lib \
                 --enable-fortran --enable-f90 \
                 --enable-omp --enable-saamg --enable-fma \
                 CC=gcc FC=gfortran F77=gfortran \
                 CFLAGS="-mcmodel=medium" CPPFLAGS="-mcmodel=medium" \
                 FCFLAGS="-mcmodel=medium" FFLAGS="-mcmodel=medium"
     make install

This has been tested under openSUSE Leap 15.2 and Linux Mint 20.1 (some
modifications might be needed under different systems).

For the **Intel compiler**, replace ‘gcc’ and ‘gfortran’ by ‘icc’ and
‘ifort’, respectively.

Installation under /opt usually requires admin rights. The same holds
for the common alternative /usr/local. For a local installation, replace
it by ‘/home/:math:`<`\ my_user_name\ :math:`>`/local’.

.. _sect_xyz:

Domain XYZ
==========

| This framework allows creating new domains (Laurentide ice sheet,
  simple testing geometry etc.). The directory
  sicopolis/src/subroutines/xyz, which hosts the domain-specific
  subroutines, is by default empty. If you want to create a new domain,
  copy the subroutines from the most similar existing domain (northern
  hemisphere, EISMINT etc.), e.g.:
| ``cp sicopolis/src/subroutines/nhem/*.F90 \``
| ``sicopolis/src/subroutines/xyz/``
| Then modify the routines according to your needs. Input files
  (topography etc.) must be placed in sicopolis/sico_in/xyz and
  specified in the run-specification header file \*.h as usual. The
  domain must be defined by the domain code ‘#define XYZ’ in the header
  file. For flexible testing, it is recommended to set the parameter
  CHECK_RES_IMAX_JMAX (compatibility check between horizontal resolution
  and number of grid points) to 0. If the new domain requires new global
  variables, they can be defined in the module
  sicopolis/src/subroutines/xyz/sico_vars.F90.

The subroutines for ISMIP HEINO are available in
sicopolis/src/subroutines/xyz/heino, and the input files are in
sicopolis/sico_in/xyz. If you copy the subroutines from
sicopolis/src/subroutines/xyz/heino to sicopolis/src/subroutines/xyz,
you can run ISMIP HEINO experiments (e.g., the run v5_heino50_st for
which a header file is available).

--------------

.. _sect_table:

Table: Simulations and computing times
======================================

.. container::
   :name: table_comp_times

   .. table:: Model times, time steps and computing (CPU) times for the
   EISMINT, Greenland and Antarctica simulations contained in the script
   multi_sico_1.sh, run with SICOPOLIS V5-dev (revision
   develop_239_rv5.2-122-g9c909c3) and the Intel Fortran Compiler 19.1
   for Linux (optimization options -xHOST -O3 -no-prec-div) on a 12-Core
   Intel Xeon Gold 6256 (3.6 GHz) PC under openSUSE Leap 15.4.
   :math:`^\dagger`: If one value is given, this is the common dynamic
   (velocity, ice thickness) and thermodynamic (temperature, water
   content, age) time step. If two values are given (marked by the
   dagger (:math:`^\dagger`) symbol), the first one is the dynamic, the
   second one the thermodynamic time step.
   :math:`^\ddagger`: All runs were done on a single core only. The
   v5_ant64_b2_xxx runs that include ice shelves can be done on multiple
   cores using OpenMP for the SSA solver. However, at the employed, low
   resolution of 64 km the solver does not scale well, and the gain in
   wall clock time by using multiple cores is very small.
   :math:`^\ast`: For this run, see the remark in
   Sect. `8.2 <#ssect_resolution_doubler>`__ on the tool
   resolution_doubler.

      +---------------+---------------+---------------+---+---------------+
      | Run           | Model time    | Time          |   |               |
      |               |               | step\ :ma     |   |               |
      |               |               | th:`^\dagger` |   |               |
      +===============+===============+===============+===+===============+
      | v5_vialov3d25 | :math:`100\   | :math:`20\    | 1 | :math:`0\;    |
      |               | ;\mathrm{ka}` | ;\mathrm{ a}` |   | \mathrm{min}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_em         | :math:`200\   | :math:`20\    | 3 | :math:`9\;    |
      | tp2sge25_expA | ;\mathrm{ka}` | ;\mathrm{ a}` |   | \mathrm{min}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_grl        | :math:`25\    | :math:`5\     | 9 | :math:`7\;    |
      | 16_bm4_ss25ka | ;\mathrm{ka}` | ;\mathrm{ a}` |   | \mathrm{min}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_an         | :math:`25\    | :math:`10\    | 5 | :math:`0\;    |
      | t40_b2_ss25ka | ;\mathrm{ka}` | ;\mathrm{ a}` |   | \mathrm{min}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_grl        | :math:`140\   | :math:`5\     | 0 | :math:`8\;    |
      | 20_b2_paleo21 | ;\mathrm{ka}` | ;\mathrm{ a}` |   | \mathrm{hrs}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_grl10      | :math:`9\     | :math:`1\     | 1 | :math:`0\;    |
      | _b2_paleo21\  | ;\mathrm{ka}` | ;\mathrm{ a}` |   | \mathrm{hrs}` |
      | :math:`^\ast` |               |               |   |               |
      +---------------+---------------+---------------+---+---------------+
      | v5_grl10_b2_  | :math:`100\   | :math:`1\     | 0 | :math:`9\;    |
      | future21_ctrl | ;\mathrm{ a}` | ;\mathrm{ a}` |   | \mathrm{min}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_grl10_b2_  | :math:`100\   | :math:`1\     | 0 | :math:`9\;    |
      | future21_asmb | ;\mathrm{ a}` | ;\mathrm{ a}` |   | \mathrm{min}` |
      +---------------+---------------+---------------+---+---------------+
      | v5_           | :math:`100\   | :math:`2\,/\, | 4 | :math:`1\;    |
      | ant64_b2_spin | ;\mathrm{ a}` | 10\;\mathrm{  |   | \mathrm{sec}` |
      | up09_init100a |               | a}^{\dagger}` |   |               |
      +---------------+---------------+---------------+---+---------------+
      | v5            | :math:`140\   | :math:`5\,/\, | 0 | :math:`7\;    |
      | _ant64_b2_spi | ;\mathrm{ka}` | 10\;\mathrm{  |   | \mathrm{hrs}` |
      | nup09_fixtopo |               | a}^{\dagger}` |   |               |
      +---------------+---------------+---------------+---+---------------+
      | v5_ant6       | :math:`500\   | :math:`2\,/\, | 0 | :math:`5\;    |
      | 4_b2_spinup09 | ;\mathrm{ a}` | 10\;\mathrm{  |   | \mathrm{min}` |
      |               |               | a}^{\dagger}` |   |               |
      +---------------+---------------+---------------+---+---------------+
      | v5_ant64_b2_  | :math:`100\   | :math:`2\,/\, | 6 | :math:`1\;    |
      | future09_ctrl | ;\mathrm{ a}` | 10\;\mathrm{  |   | \mathrm{sec}` |
      |               |               | a}^{\dagger}` |   |               |
      +---------------+---------------+---------------+---+---------------+

--------------
