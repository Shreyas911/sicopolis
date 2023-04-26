.. _domain_grid_time:

Model domain, grid and time
***************************

.. _model_domain:

Model domain
============

.. _defined_domain:

Selecting a pre-defined domain
------------------------------

SICOPOLIS provides several pre-defined model domains. They can be chosen by defining a domain code in the run-specs header as follows\:

.. code-block:: fortran

  #define GRL
  !             Simulated domain:
  !               ANT     - Antarctica
  !               ASF     - Austfonna
  !               EISMINT - EISMINT (Phase 2 SGE and modifications)
  !               GRL     - Greenland
  !               NHEM    - Northern hemisphere
  !               SCAND   - Scandinavia
  !               TIBET   - Tibet
  !               NMARS   - North polar cap of Mars
  !               SMARS   - South polar cap of Mars

This example would select the domain for the Greenland ice sheet. Correspondingly for the other listed domains.

.. _new_domain:

Setting up a new domain
-----------------------

In addition to the pre-defined domains, there is an unspecified domain XYZ. This framework allows creating new domains (Laurentide ice sheet, some ice cap, simple testing geometry...) quite easily. The directory src/subroutines/xyz, which hosts the domain-specific subroutines, is by default empty. If you want to create a new domain, copy the subroutines from the most similar existing domain, e.g., starting from Antarctica::

  cp src/subroutines/ant/*.F90 src/subroutines/xyz/

Then modify the routines according to your needs. Input files (topography etc.) must be placed in sico_in/xyz and specified in the run-specs header file \*.h as usual. The domain must be defined by the domain code

.. code-block:: fortran

  #define XYZ

in the header. For flexible testing, it is recommended to deactivate the compatibility check between horizontal resolution and number of grid points\:

.. code-block:: fortran

  #define CHECK_RES_IMAX_JMAX 0

If the new domain requires new global variables, they can be defined in the module src/subroutines/xyz/sico_vars.F90.

The subroutines for ISMIP HEINO (Calov et al. :cite:`calov_etal_2010`) are available in src/subroutines/xyz/heino, and the input files are in sico_in/xyz. If you copy the subroutines from src/subroutines/xyz/heino to src/subroutines/xyz, you can run ISMIP HEINO experiments (e.g., the run v5_heino50_st for which a run-specs header file is available).

.. _spatial_grid:

Spatial grid
============

In principle, SICOPOLIS allows using any orthogonal coordinates on the Earth's surface, provided that the two components :math:`g_{11}` and :math:`g_{22}` of the metric tensor are known (see ":ref:`orthog_coord`"). They are computed in the module ``metric_m``. Three options are currently implemented and can be selected in the run-specs header\:

.. code-block:: fortran

  #define GRID 1
  !                0 : Cartesian coordinates in the stereographic plane
  !                    without distortion correction
  !                1 : Cartesian coordinates in the stereographic plane
  !                    with distortion correction
  !                2 : Geographical coordinates (longitude/latitude)

For the most common case of Cartesian coordinates :math:`x` and :math:`y` in the stereographic plane (or any other projection plane), let the domain be the rectangle described by :math:`[x_0,x_\mathrm{max}]`, :math:`[y_0,y_\mathrm{max}]`. It is discretized by a regular (structured) grid with horizontal resolution :math:`\Delta{x}`, which is the same for the :math:`x`- and :math:`y`-directions. The location of the grid points :math:`x_i` and :math:`y_j` is then given by

.. math::
  :label: eq_discr_x

  x_i = x_0 + i\Delta{x}, \qquad i=0\,(1)\,i_\mathrm{max},

.. math::
  :label: eq_discr_y

  y_j = y_0 + j\Delta{x}, \qquad j=0\,(1)\,j_\mathrm{max},

where the notation :math:`a\,(b)\,c` means "from :math:`a` to :math:`c` in steps of :math:`b`". Note that the indices :math:`i` and :math:`j` run from 0, so that the number of grid points is actually :math:`i_\mathrm{max}+1` and :math:`j_\mathrm{max}+1`, respectively. In the run-specs headers, the parameters to be defined are

* X0 (:math:`=x_0`, :math:`x` coordinate of the origin point in km),
* Y0 (:math:`=y_0`, :math:`y` coordinate of the origin point in km),
* DX (:math:`=\Delta{}x`, horizontal grid spacing in km),
* IMAX (:math:`=i_\mathrm{max}`, maximum value of the index :math:`i`),
* JMAX (:math:`=j_\mathrm{max}`, maximum value of the index :math:`j`).

For the vertical (:math:`z`) direction, a terrain-following ("sigma") transformation is employed that maps vertical columns in the physical space onto :math:`[0,1]` intervals. If the polythermal two-layer method (POLY, see Section ":ref:`ice_thermodynamics`") is employed, this mapping is done separately for the upper cold-ice layer (:math:`\zeta_\mathrm{c}` domain), the lower temperate-ice layer (:math:`\zeta_\mathrm{t}` domain) and the lithosphere layer (:math:`\zeta_\mathrm{r}` domain). The transformation is linear for the :math:`\zeta_\mathrm{t}` and :math:`\zeta_\mathrm{r}` domains. However, for the :math:`\zeta_\mathrm{c}` domain, exponential stretching is used so that equidistant grid points in the transformed domain map on grid points concentrating towards the base in the physical :math:`z`-coordinate\:

.. math::
  :label: eq_sigma_trans

  \mbox{Transformation equations for }
  \zeta_\mathrm{c},\; \zeta_\mathrm{t},\; \zeta_\mathrm{c}\; ...

The location of the grid points in the three transformed domains is given by

.. math::
  :label: eq_discr_zc

  (\zeta_\mathrm{c})_{k_\mathrm{c}} = k_\mathrm{c}/k_\mathrm{c,max},
  \qquad k_\mathrm{c}=0\,(1)\,k_\mathrm{c,max},

.. math::
  :label: eq_discr_zt

  (\zeta_\mathrm{t})_{k_\mathrm{t}} = k_\mathrm{t}/k_\mathrm{t,max},
  \qquad k_\mathrm{t}=0\,(1)\,k_\mathrm{t,max},

.. math::
  :label: eq_discr_zr

  (\zeta_\mathrm{r})_{k_\mathrm{r}} = k_\mathrm{r}/k_\mathrm{r,max},
  \qquad k_\mathrm{r}=0\,(1)\,k_\mathrm{r,max}.

The number of grid points results as :math:`k_\mathrm{c,max}+1`, :math:`k_\mathrm{t,max}+1` and :math:`k_\mathrm{r,max}+1`. The parameters in the run-specs headers are

* KCMAX (:math:`=k_\mathrm{c,max}`, maximum value of the index :math:`k_\mathrm{c}`),
* KTMAX (:math:`=k_\mathrm{t,max}`, maximum value of the index :math:`k_\mathrm{t}`),
* KRMAX (:math:`=k_\mathrm{r,max}`, maximum value of the index :math:`k_\mathrm{r}`),
* DEFORM (:math:`=a`, exponential stretch parameter for the :math:`k_\mathrm{c}` domain; 0.0d0 produces a linear transformation).

For all other thermodynamics schemes (COLD, ENTC, ENTM; see Section ":ref:`ice_thermodynamics`"), ...

Topography...

.. _model_time:

Model time
==========

Initial time, final time, time steps...

.. math::
  :label: eq_discr_t

  t^n = t^0 + n\Delta{}t, \qquad n=0\,(1)\,n_\mathrm{max}.

Lorem ipsum...
