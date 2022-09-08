.. _tapenade:

Open source AD Tool TAPENADE
****************************

TAPENADE is an Automatic Differentiation Engine developed at Inria at Sophia Antipolis by the Tropics then Ecuador teams. TAPENADE takes as input a computer source program,
plus a request for differentiation. TAPENADE builds and returns the differentiated source program, that evaluates the required derivatives.

.. _build_tapenade:

Building TAPENADE (latest 3.16)
===============================

While the SICOPOLIS source files are prepared to generate adjoint sensitivities, they will not be able to do so without an operable installation of Tapenade. Fortunately the Tapenade installation procedure is straight forward.

We detail the instructions here, but the latest instructions can always be found `here. <http://www-sop.inria.fr/ecuador/tapenade/distrib/README.html>`__

Prerequisites for Linux or Mac OS X
-----------------------------------

Before installing Tapenade, you must check that an up-to-date Java Runtime Environment is installed. Tapenade will not run with older Java Runtime Environment.

Steps for Linux or Mac OS X
---------------------------

1. Read `the Tapenade license. <https://tapenade.gitlabpages.inria.fr/userdoc/build/html/LICENSE.html>`__

2. Download `tapenade_3.16.tar <http://www-sop.inria.fr/ecuador/tapenade/distrib/tapenade_3.16.tar>`__ into your chosen installation directory *install_dir*.

**NOTE**: Alternatively, a Tapenade version that works correctly with SICOPOLIS-AD v2 is always available in the ``test_ad/tapenade_3.16`` directory.

3. Go to your chosen installation directory *install_dir*, and extract Tapenade from the tar file :

::

    % tar xvfz tapenade_3.16.tar

4. On Linux, depending on your distribution, Tapenade may require you to set the shell variable ``JAVA_HOME`` to your java installation directory. It is often ``JAVA_HOME=/usr/java/default``. You might also need to modify the ``PATH`` by adding the bin directory from the Tapenade installation. An example can be found :ref:`here <tapenade_bashrc_snippet>`.

Prerequisites for Windows
-------------------------

**NOTE**: Although Tapenade can be built on Windows, SICOPOLIS requires a Unix-like system (e.g., Linux), as mentioned :ref:`here <sico_prerequisites>`.

Before installing Tapenade, you must check that an up-to-date Java Runtime Environment is installed. Tapenade will not run with older Java Runtime Environment. The Fortran parser of Tapenade uses `cygwin <https://www.cygwin.com/>`__.

Steps for Windows
-----------------

1. Read `the Tapenade license. <https://tapenade.gitlabpages.inria.fr/userdoc/build/html/LICENSE.html>`__

2. Download `tapenade_3.16.zip <http://www-sop.inria.fr/ecuador/tapenade/distrib/tapenade_3.16.zip>`__ into your chosen installation directory *install_dir*.

**NOTE**: Alternatively, a Tapenade version that works correctly with SICOPOLIS-AD v2 is always available in the ``test_ad/tapenade_3.16`` directory.

3. Go to your chosen installation directory *install_dir*, and extract Tapenade from the zip file.

4. Save a copy of the ``install_dir\tapenade_3.16\bin\tapenade.bat`` file and modify ``install_dir\tapenade_3.16\bin\tapenade.bat`` according to your installation parameters:

replace ``TAPENADE_HOME=..`` by ``TAPENADE_HOME="install_dir"\tapenade_3.16``
replace ``JAVA_HOME="C:\Progra~1\Java\jdkXXXX"`` by your current java directory
replace ``BROWSER="C:\Program Files\Internet Explorer\iexplore.exe"`` by your current browser.

.. _tapenade_bashrc_snippet:

**NOTE**: Every time you wish to use the adjoint capability of SICOPOLIS-AD, you must re-source the environment. We recommend that this be done automatically in your bash or c-shell profile upon login. An example of an addition to a ``.bashrc`` file from a Linux server is given below. Luckily, shell variable ``JAVA_HOME`` was not required to be explicitly set for this particular Linux distribution, but might be necessary for some other distributions.

::

    ##set some env variables for SICOPOLIS tapenade

    export TAPENADE_HOME="/home/shreyas/tapenade_3.16"
    export PATH="$PATH:$TAPENADE_HOME/bin"

    ##Modules

    module use /share/modulefiles/
    module load java/jdk/16.0.1 # Java required by Tapenade


You should now have a working copy of Tapenade.

For more information on the tapenade command and its arguments, type :

::

    tapenade -?

