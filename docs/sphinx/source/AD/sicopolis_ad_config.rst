.. _sicopolis_ad_config:

Ice sheet model SICOPOLIS
*************************

.. _sico_prerequisites:

Prerequisites for SICOPOLIS
===========================

See the ":ref:`getting_started`" section. However, we will mention here some steps since using the Automatic Differentiation (AD) capabilities with Tapenade requires a slightly modified setup.

The satisfaction of the following prerequisites is highly recommended to access all the features of the code. Details can differ from the ":ref:`getting_started`" section, since there are multiple ways to do things. We detail one of them here.

GNU GCC Compiler (gfortran+gcc) or Intel Compiler (ifort+icc)
-------------------------------------------------------------

We have tested the software on gfortran/gcc v5.4.0, v7.2.0 and v8.5.0, any intermediate versions should work just as well. We have also tested the software on ifort/icc v18.0.0 (however, it should be noted that we have not tested the external Lis solver with Intel compilers).

Git
---

The Git repository of SICOPOLIS is kindly hosted by the GitLab system of the Alfred Wegener Institute for Polar and Marine Research (AWI) in Bremerhaven, Germany. Front page `here <https://gitlab.awi.de/sicopolis/sicopolis/>`__.

Cloning the latest develop or ad revision::

  git clone --branch develop \
  https://gitlab.awi.de/sicopolis/sicopolis.git

  git clone --branch ad \
  https://gitlab.awi.de/sicopolis/sicopolis.git

(Cloning with SSH instead of HTTPS is also available. See the above GitLab front page link for details.)

You should then have a new directory ``sicopolis`` that contains the entire program package.

Lis (1.4.43 or newer)
---------------------

Lis installation is mandatory to use shallow-shelf/shelfy-stream dynamics in simulations. Install Lis as explained in :ref:`Dependencies/Lis <dependencies-lis>`. The following commands might be helpful, they are written for the latest version at the time of writing::

  wget https://www.ssisc.org/lis/dl/lis-2.0.30.zip
  unzip lis-2.0.30.zip
  cd lis-2.0.30
  ./configure --enable-fortran --prefix=/lis-2.0.30/installation/ --enable-shared && make && make check && make install

For AD purposes, we compile the code using the ``src/MakefileTapenade`` makefile. This makefile requires the following environment variables set\:

1. ``LISDIR`` - The installation directory for the Lis version to be used. You can change this variable, either in the Makefile directly, or automatically in your bash or c-shell profile upon login (for example, ``.bashrc``). Examples for both are shown here.

  ``src/MakefileTapenade``::

    ifndef LISDIR
    LISDIR=/home/shreyas/lis-2.0.30/installation
    endif

  ``.bashrc``::

    export LISDIR="/home/shreyas/lis-2.0.30/installation"

2. ``LIBLIS`` - Absolute path to ``liblis.so``. By default in ``src/MakefileTapenade``, it is ``${LISDIR}/lib/liblis.so``. If you follow the original instructions to install Lis, this should work, else one can set it manually within ``src/MakefileTapenade``. 

3. ``LIBLISFLAG`` - Add directory using ``-L`` to be searched for ``-llis``. By default in ``src/MakefileTapenade``, it is ``-L${LISDIR}/lib -llis``. If you follow the original instructions to install Lis, this should work, else one can set it manually within ``src/MakefileTapenade``.

4. ``LISFLAG`` - This flag declares directories to be searched for Lis ``#include`` header file ``lisf.h``, as well as defines the ``BUILD_LIS`` as a macro with value 1. By default in ``src/MakefileTapenade``, it is ``-DBUILD_LIS -I${LISDIR}/include/``. If you follow the original instructions to install Lis, this should work, else one can set it manually within ``src/MakefileTapenade``.

**NOTE**: Some users have reported needing to extend their ``LD_LIBRARY_PATH`` with the location of ``${LISDIR}/lib`` in order to find ``liblis.so.0``.

NetCDF (3.6.x or newer)
-----------------------

NetCDF installation is mandatory since it is a powerful library with widespread use for I/O with a machine-independent data format. Install NetCDF as explained in :ref:`Dependencies/NetCDF <dependencies-netcdf>`. In some cases, for example while working on a shared server which uses a module manager or Docker container, thing have to be set up differently. ``src/MakefileTapenade`` needs either the ``NETCDF_FORTRAN_DIR`` macro set or both ``NETCDF_F90_FLAG`` and ``LIB_NETCDF_F90_FLAG`` set (see code snippet from ``src/MakefileTapenade`` here). ::

  ifndef NETCDF_F90_FLAG
  ifndef LIB_NETCDF_F90_FLAG
  ifdef NETCDF_FORTRAN_DIR
  LIB_NETCDF_F90=${NETCDF_FORTRAN_DIR}/lib/libnetcdff.so
  LIB_NETCDF_F90_FLAG=-L${NETCDF_FORTRAN_DIR}/lib -lnetcdff
  NETCDF_F90_FLAG=-I${NETCDF_FORTRAN_DIR}/include/
  endif
  endif
  endif

1. ``NETCDF_FORTRAN_DIR`` - The installation directory for netcdf-fortran. You can change this variable, either in the Makefile directly, or automatically in your bash or c-shell profile upon login (for example, ``.bashrc``). Examples for both are shown here.

  ``src/MakefileTapenade`` ::

    ifndef NETCDF_FORTRAN_DIR
    NETCDF_FORTRAN_DIR=/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4
    endif

  ``.bashrc`` ::

    export NETCDF_FORTRAN_DIR="/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4"

2. ``LIB_NETCDF_F90`` - Absolute path to ``libnetcdff.so``. By default in ``src/MakefileTapenade``, it is ``{NETCDF_FORTRAN_DIR}/lib/libnetcdff.so``.

3. ``LIB_NETCDF_F90_FLAG`` - Add directory using ``-L`` to be searched for ``-lnetcdff``. By default in ``src/MakefileTapenade``, it is ``-L${NETCDF_FORTRAN_DIR}/lib -lnetcdff``. See some examples below where this has to be set explicitly in case of a Docker container.

4. ``NETCDF_F90_FLAG`` - This flag declares directories to be searched for netcdf-fortran ``#include`` header files. By default in ``src/MakefileTapenade``, it is ``-I${NETCDF_FORTRAN_DIR}/include/``. See some examples below where this has to be set explicitly in case of a Docker container. 

For a server that uses modules, you can load the relevant modules using commands like these (can also make permanent by adding to login script like ``.bashrc``::

  % module use /share/modulefiles/
  % module load openmpi
  % module load netcdf-fortran
  % module load netcdf

You then have to give the ``NETCDF_FORTRAN_DIR`` macro to ``src/MakefileTapenade``, either by adding to a login script or directly inside the makefile. If your system uses a module manager, you can query to find the exact directory location. ::

  % module show netcdf-fortran
  ----------------------------------------------------------
  /opt/ohpc/pub/moduledeps/gnu-openmpi/netcdf-fortran/4.4.4:
  ----------------------------------------------------------

  whatis("Name: NETCDF_FORTRAN built with gnu toolchain ")
  whatis("Version: 4.4.4 ")
  whatis("Category: runtime library ")
  whatis("Description: Fortran Libraries for the Unidata network Common Data Form ")
  whatis("http://www.unidata.ucar.edu/software/netcdf/ ")
  prepend_path("PATH","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/bin")
  prepend_path("MANPATH","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/share/man")
  prepend_path("INCLUDE","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/include")
  prepend_path("LD_LIBRARY_PATH","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/lib")
  setenv("NETCDF_FORTRAN_DIR","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4")
  setenv("NETCDF_FORTRAN_BIN","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/bin")
  setenv("NETCDF_FORTRAN_LIB","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/lib")
  setenv("NETCDF_FORTRAN_INC","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/include")
  help([[ 
  This module loads the NetCDF Fortran API built with the gnu compiler toolchain.
   
  Note that this build of NetCDF leverages the HDF I/O library and requires linkage
  against hdf5 and the native C NetCDF library. Consequently, phdf5 and the standard C
  version of NetCDF are loaded automatically via this module. A typical compilation
  example for Fortran applications requiring NetCDF is as follows:
   
  ]])

In this case ``NETCDF_FORTRAN_DIR=/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4``.

For a Docker container, for example one with a ``centos:8`` distribution, and the ``dnf`` package manager, NetCDF is typically installed as follows::

  RUN dnf install -y https://github.com/openhpc/ohpc/releases/download/v1.3.GA/ohpc-release-1.3-1.el7.x86_64.rpm
  
  # Add some packages
  RUN dnf -y install epel-release
  RUN dnf -y install dnf-plugins-core
  RUN dnf config-manager --set-enabled powertools
  RUN dnf -y install make which git
  RUN dnf -y install diffutils
  RUN dnf -y install vim
  RUN dnf -y install autoconf automake
  RUN dnf -y install valgrind-ohpc
  RUN dnf -y install gnu8-compilers-ohpc
  RUN dnf -y install gsl-gnu8-ohpc hdf5-gnu8-ohpc
  RUN dnf -y install openmpi-devel
  RUN dnf -y install bc wget zlib-devel perl-Digest-MD5
  RUN dnf -y --enablerepo=powertools install netcdf-fortran netcdf-devel # NetCDF installation
  RUN dnf -y install netcdf-fortran-devel # NetCDF installation

In this case, you will find that the ``./usr/lib64/gfortran/modules/netcdf.mod`` exists in your docker environment. In this case, you can just directly set  ``NETCDF_F90_FLAG=-I/usr/lib64/gfortran/modules`` either the makefile or the login script (no need to set ``NETCDF_FORTRAN_DIR`` macro). 

You can also confirm that the files ``/usr/lib64/libnetcdff.so*`` and ``/usr/lib64/libnetcdf.so*`` exist, which means you have to set ``LIB_NETCDF_F90_FLAG=-L/usr/lib64 -lnetcdff``.

The instructions given in :ref:`Dependencies/NetCDF <dependencies-netcdf>`, and these two cases should help cover most of the issues with the installation of NetCDF.

Unix-like system
----------------

A Unix-like system, e.g. Linux (Ubuntu, CentOS, Fedora, Redhat, etc.), MacOS is required to run both SICOPOLIS and SICOPOLIS-AD v2.

Setting up SICOPOLIS
====================

The Git repository of SICOPOLIS is kindly hosted by the GitLab system of the `Alfred Wegener Institute for Polar and Marine Research (AWI) <https://www.awi.de/>`__ in Bremerhaven, Germany. 

* Front page: `Front page: https://gitlab.awi.de/sicopolis/sicopolis/ <https://gitlab.awi.de/sicopolis/sicopolis/>`__

* Cloning the latest ``ad`` (the branch most relevant to us) revision with Git::

    git clone --branch ad \
    https://gitlab.awi.de/sicopolis/sicopolis.git

  Cloning with SSH instead of HTTPS is also available. See the above GitLab link for details.

* Tagged versions of SICOPOLIS can be accessed from the `archive <http://www.sicopolis.net/archive/>`__.

SICOPOLIS and SICOPOLIS-AD v2 applications are built using a configuration header file in ``runs/headers``. A typical user setup involves copying over example configuration files from ``runs/headers/templates`` (see below), and suitably modifying one of them for custom runs.

Initial configuration
===================== 

In addition to the steps above, the following steps need to be performed from the root of the repository\:

* Copy template header files from ``runs/headers/templates`` to ``runs/headers`` so that SICOPOLIS can read one of these header files for the simulations desired by the user. Also, one can modify them suitably for their own custom simulations. The original files are always stored in ``runs/headers/templates`` for reference. Run the following command from the root directory of the repository::

    ./copy_templates.sh

* Get the input data files needed for both Greenland and Antarctica. These files are stored on a server and needed for various inputs such as geothermal heat flux, physical parameters, height of the ice base and lithosphere, precipitation, definition of regions for heterogenous basal sliding, etc. Run the following command from the root directory of the repository::

    ./get_input_files.sh

* Locate the file ``sico_environment.sh`` in the directory ``sicopolis/runs``, open it with a text editor, and replace the "Default" entry for ``SICO_INSTITUTION`` by the name of your institution (max. 256 characters). This is just for bookkeping purposes.

Now, you are ready to use SICOPOLIS-AD v2, as described in :ref:`Running SICOPOLIS-AD v2 <running>`!
