[![Documentation Status](https://readthedocs.org/projects/sicopolis/badge/?version=latest)](https://sicopolis.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://gitlab.awi.de/sicopolis/sicopolis/badges/ad/pipeline.svg)](https://gitlab.awi.de/sicopolis/sicopolis/-/commits/ad)

SICOPOLIS-AD v2
---------------

`SICOPOLIS-AD v2` is an open-source adjoint and tangent linear modeling framework for ice sheet model `SICOPOLIS` enabled by the Automatic Differentiation (AD)  tool `Tapenade`.

Find the documentation [here](https://sicopolis.readthedocs.io/en/latest/).
Detailed instructions for installation are given [here](https://sicopolis.readthedocs.io/en/latest/AD/AD.html#configuration)

The AD code is tested using `Gitlab-CI`. The testing setup can be found [here](test_ad/). 

`SICOPOLIS-AD v2` is an open source project that relies on the participation of its users, and we welcome contributions. Users can contribute using the usual pull request mechanisms in git, and if the contribution is substantial, they can contact us to discuss gaining direct access to the repository.

If you think you’ve found a bug, please check if you’re using the latest version of the model. If the bug is still present, then think about how you might fix it and file a ticket in the Gitlab issue tracker (you might need to request membership access on Gitlab, which we can approve). Your ticket should include: what the bug does, the location of the bug: file name and line number(s), and any suggestions you have for how it might be fixed.

To request a new feature, or guidance on how to implement it yourself, please open a ticket with a clear explanation of what the feature will do.

You can also directly contact Shreyas Gaikwad (shreyas.gaikwad@utexas.edu) for any of the above.

Installation
------------

* Tapenade 3.16 (or latest version)

* SICOPOLIS v5.3

	- Pre-requisites for smooth functioning

		1. GNU GCC compiler (gfortran, tested on 5.4.0 or newer)

		2. Git (for cloning the SICOPOLIS repo)

		3. Unix-like system

		4. LIS (External library similar to PETSc, 1.4.43 or newer)

		5. NetCDF (External library and machine-independent data format, 3.6.x or newer)

For more details, check the documentation [here](https://sicopolis.readthedocs.io/en/latest/AD/AD.html#installation). 

Capabilities
------------

* Gradient calculation enabled using
	1. Adjoint mode
	2. Tangent linear mode
	3. Finite differences mode 

* Automated python scripts for adjoint validation with finite differences

Potential Uses
--------------

* Sensitivity analysis

* Model Calibration

* Efficient uncertainty quantification

* Optimal Experimental Design (OED)

Contributors
------------

Shreyas Gaikwad, Sri Hari Krishna Narayanan, Laurent Hascoet, Liz Curry-Logan, Ralf Greve, Patrick Heimbach

For conributions of each author, please check the documentation [here](https://sicopolis.readthedocs.io/en/latest/AD/AD.html#contributors).

License
-------

The software is hosted under the GNU General Public License.
