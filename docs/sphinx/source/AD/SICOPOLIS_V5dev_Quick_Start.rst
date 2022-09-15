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
