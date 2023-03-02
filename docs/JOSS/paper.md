---
title: 'SICOPOLIS-AD v2: tangent linear and adjoint modeling framework for ice sheet modeling enabled by automatic differentiation tool Tapenade'

tags:
  - SICOPOLIS-AD
  - Fortran
  - Automatic Differentiation
  - Inverse modeling
  - Data assimilation

authors:
  - name: Shreyas Sunil Gaikwad^[Corresponding author]
    orcid: 0000-0003-2079-4218
    affiliation: 1
  - name: Laurent Hascoet
    orcid: 0000-0002-5361-0713
    affiliation: 2
  - name: Sri Hari Krishna Narayanan
    orcid: 0000-0003-0388-5943
    affiliation: 3
  - name: Liz Curry-Logan
    affiliation: 1
  - name: Ralf Greve
    orcid: 0000-0002-1341-4777
    affiliation: "4, 5" # (Multiple affiliations must be quoted)
  - name: Patrick Heimbach
    orcid: 0000-0003-3925-6161
    affiliation: "1, 6, 7" # (Multiple affiliations must be quoted)

affiliations:
 - name: Oden Institute for Computational Engineering and Sciences, University of Texas at Austin, USA
   index: 1
 - name: Institut National de Recherche en Informatique et Automatique, France
   index: 2
 - name: Mathematics and Computer Science Division, Argonne National Laboratory, USA
   index: 3
 - name: Institute of Low Temperature Science, Hokkaido University, Japan
   index: 4
 - name: Arctic Research Center, Hokkaido University, Japan
   index: 5
 - name: Jackson School of Geosciences, University of Texas at Austin, USA
   index: 6
 - name: Institute for Geophysics, University of Texas at Austin, USA
   index: 7

date:
bibliography: paper.bib
---

# Summary

SImulation COde for POLythermal Ice Sheets ([`SICOPOLIS`](https://gitlab.awi.de/sicopolis/sicopolis/-/tree/master)) is an open-source, 3D dynamic/thermodynamic model that simulates the evolution of large ice sheets and ice caps. SICOPOLIS has been developed continuously and applied to problems of past, present, and future glaciation of Greenland, Antarctica, and others. It uses the finite differences discretization on a staggered Arakawa C grid and employs the shallow ice and shallow shelf approximations, making it suitable for paleoclimatic simulations. We present a new framework for generating derivative code, i.e., tangent linear, adjoint, or Hessian models, of SICOPOLIS. These derivative operators are powerful computational engines to efficiently compute comprehensive gradients or sensitivities of scalar-valued model output, including least-squares model-data misfits or important quantities of interest, to high-dimensional model inputs (such as model initial conditions, parameter fields, or boundary conditions). The new version 2 ([`SICOPOLIS-AD v2`](https://gitlab.awi.de/sicopolis/sicopolis/-/tree/ad)) framework is based on the source-to-source automatic differentiation (AD) tool [`Tapenade`](https://tapenade.gitlabpages.inria.fr/userdoc/build/html/tapenade/faq.html) which has recently been open-sourced. The switch from a previous AD tool (`OpenAD`) used in SICOPOLIS-AD version 1 to `Tapenade` overcomes several limitations outlined here. The framework is integrated with the SICOPOLIS model's main trunk and is freely available.  

# Statement of need

The two contemporary ice sheets, Greenland and Antarctica, are dynamic entities whose evolution is governed by a set of nonlinear partial differential equations (PDEs) that describe the conservation of mass, momentum, and energy, as well as constitutive laws for the material properties of ice. In general, these equations cannot be solved analytically but must be solved numerically. Ice sheet models are a computer representation of these PDEs. They require as input parameters (i) initial conditions of the state of the ice sheet, (ii) surface boundary conditions, such as precipitation, (iii) basal boundary conditions, such as geothermal flux, and (iv) model parameters, such as flow law parameters.  Despite advances in numerical modeling of ice sheets, the effects of ad-hoc initialization and the uncertainties in these independent input parameters propagate to quantities of interest (QoI), such as future projections of sea-level rise, which is of economic and societal importance [@Schinko_2020]. It is thus desirable to evaluate the sensitivities of our QoI to these independent input variables.

In the context of ice sheet modeling, sensitivities of model-data misfits or other QoI are a key ingredient for performing model calibration, state estimation, or uncertainty quantification (UQ), which guide the improvement of model simulations through PDE-constrained gradient-based optimization. 

`SICOPOLIS-AD v2` leverages the recently open-sourced AD tool `Tapenade` [@hascoetTAP] to generate code for the adjoint model of the open-source ice sheet model, `SICOPOLIS` [@greve1997application; @greve2009dynamics; @greve2011initial]. Sensitivities can be calculated using a single forward and adjoint model evaluation, instead of the $\mathcal{O}(N)$ forward model evaluations. Empirically, one adjoint model evaluation is about 5-10 times as expensive as a forward model run. The adjoint computation is highly efficient for calculating sensitivities when $N$ is large (typically, $N \sim 10^4-10^6$).

The functionality to generate a tangent linear version of the forward model is also included, which was not available in `SICOPOLIS-AD v1`. This is valuable for UQ of the inferred parameters, as well as uncertainty propagation to QoIs. It can also be used to verify the results of the adjoint model. 

# Target Audience

This package is intended as a resource that enables sensitivity analysis, model calibration, and uncertainty quantification of a continent-scale ice sheet model. Our package is also intended to serve as a guide for future work in the application of open-source AD tools for physics-based simulation codes written in `Fortran`. 

# State of the field

`SICOPOLIS` is among the early thermo-mechanical models to simulate contemporary and paleo continental-scale ice sheets [@greve1997application]. Like similar models developed at the time, including Glimmer and its successor, the Community Ice Sheet Model (`CISM`) [@https://doi.org/10.1029/2008JF001015], `GRISLI` [@Ritz1996], the model by @huybrechts19903, or by @pollard2009, `SICOPOLIS` has been based (until recently) on the so-called shallow ice approximation to simplify the Cauchy stress tensor in the momentum conservation equation, implemented on a regular, finite-difference mesh. See @https://doi.org/10.1029/2003JF000065 for other approximations commonly employed in ice sheet models. This approximation enabled the efficient computation of ice sheet evolution over long, glacial/deglacial cycles.

The last decade has seen substantial advances in continental-scale ice sheet modeling, with the development of several new ice sheet models (some of which are on unstructured grids using finite element methods), notably the Ice Sheet System Model `ISSM` [@https://doi.org/10.1029/2011JF002140], the Parallel Ice Sheet Model `PISM` [@bueler_brown_lingle_2007], Elmer/Ice [@gmd-6-1299-2013], or the MPAS-Albany Land Ice `MALI` [@gmd-11-3747-2018]. While designed to capture the evolution of short-term, fast-flowing, or fast-changing outlet glaciers via horizontal stress contributions, these models have so far found little application in paleo-ice sheet simulations due to their extensive computational costs. A compilation of the suite of ice sheet models used for the latest Ice Sheet Model Intercomparison Project, Phase 6 (`ISMIP6`) in support of the IPCC's Sixth Assessment Report is available in @https://doi.org/10.1029/2020GL091741 and @gmd-9-4521-2016.

Relevant to this paper, of all the time-evolving models listed, apart from `SICOPOLIS-AD` [@heimbach2009greenland; @gmd-13-1845-2020], only the `ISSM` model and variants thereof possess adjoint model codes which have been generated, in part, using automatic differentiation [@tc-8-2335-2014; @doi:10.1080/10556788.2017.1396600]. Multi-centennial and longer integrations with the adjoint model have so far been conducted only with `SICOPOLIS-AD`.

# Features 

AD tools such as the commercial `TAF` [@giering_kaminski] and the open-source `OpenAD` [@openad] have been used previously with `SICOPOLIS` [@heimbach2009greenland; @gmd-13-1845-2020]. `OpenAD` is no longer actively developed because it is based on the `Open64` compiler which ceased development in 2011. The differentiation of `SICOPOLIS`, therefore, must be performed using a different tool. Compared to `OpenAD`, the `Tapenade` enabled implementation has the following advantages:

* It is up-to-date with the latest `SICOPOLIS` code.

* The AD tool `Tapenade` is open-source and actively maintained.

* A new tangent linear code generation capability is introduced.

* The AD-generated codes can accept [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) inputs.

* The external library [LIS](https://www.ssisc.org/lis/index.en.html), its tangent linear code, and adjoint code are correctly incorporated which can improve the simulation of Antarctic ice shelves and Greenland outlet glaciers.

* [Gitlab-CI](https://docs.gitlab.com/ee/ci/), a [Docker](https://www.docker.com/), and the [pytest](https://docs.pytest.org/en/7.1.x/) framework are leveraged for Continuous Integration (CI) to track changes in the trunk that "break" the AD-based code generation.

* The entire code is parsed by `Tapenade`, preventing cumbersome manual maintenance of subroutines to initialize the adjoint runs.

* Python scripts are provided for quick setup of the compilation, I/O, and execution processes based on user-provided metadata.

* The setup is [well-documented](https://sicopolis.readthedocs.io/en/latest/AD/AD.html), along with tutorials.

# Software requirements and external usage

`SICOPOLIS-AD v2` is built on top of the ice sheet model `SICOPOLIS` and uses `Tapenade` to differentiate this model. All the prerequisites of using `SICOPOLIS` and `Tapenade` need to be satisfied. A Python installation is needed to use the automation tools.

# Example 

We illustrate the use of our tool with the example of a steady-state simulation of the Greenland ice sheet under modern climate conditions. The corresponding SICOPOLIS configuration header file, `v5_grl16_bm5_ss25ka`, is provided as a reference template in the standard `SICOPOLIS` distribution. We shorten the total integration time to 100 simulated years to keep the computational cost of the tangent linear and finite differences reasonable. Our QoI (i.e., dependent variable) is the total volume of the ice sheet at the end of the run (`fc`). The sensitivity is evaluated with respect to the geothermal heat flux, `q_geo` (independent variable), a 19,186-dimensional field. The results are shown in Figure 1. 

![Validation exercise for adjoint (ADM) and tangent linear (TLM) models using the finite differences (FD) results for the sensitivity of `fc` with respect to `q_geo`. The upper row shows the sensitivities computed using the adjoint model (reverse-mode AD), tangent linear (forward-mode AD), and finite differences, respectively. The bottom run illustrates the relative error between (ADM, FD), (TLM, FD), and (ADM, TLM) respectively. For the bottom row, note that the values of relative error are only shown for points where the value of the gradient is "significant", i.e. within 4 orders of magnitude of the maximum absolute value of the gradient. \label{fig:fig1}](Figure_1.png)

The results show good agreement between all three modes used to evaluate this sensitivity. The error is less than $6\%$ between AD-generated (adjoint/tangent linear codes) and finite differences at all but one point with "significant" gradient values, i.e. within 4 orders of magnitude of the maximum absolute value of the finite differences gradient. The relative error between the AD-generated adjoint and tangent linear models is less than $0.002\%$ at all points with values within 4 orders of magnitude of the maximum absolute value of the finite differences gradient. However, the adjoint model is much faster than the other two, as shown in Table 1, because the number of evaluations of the latter two scales linearly with the parameter dimension (~$\mathcal{O}(N)$). The discrepancy will be even larger if a finer mesh is used.

| Gradient calculation method | Time (in seconds) for 16 km mesh |
| :-------------------------: | :------------------------------: |
| Finite Differences | $1.640 \times 10^5$ |
| Tangent Linear Model | $9.793 \times 10^4$ |
| Adjoint Model | $2.214 \times 10^1$ | 

: Comparison of the time taken by various methods to evaluate the gradient for a scalar objective function with respect to a 19,186-dimensional 2D field (16 km mesh) in a typical SICOPOLIS run. The runs are performed on Intel Xeon CPU E5-2695 v3 nodes (2.30 GHz clock rate, 35.84 MB L3 cache, 63.3 GB memory).

# Acknowledgements

This work was supported by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research Program, under contract number DE-AC02-06CH11357, by National Science Foundation OPP/P2C2 grant #1903596, by Japan Society for the Promotion of Science (JSPS) KAKENHI grant numbers JP17H06104 and JP17H06323, and by the Japanese Ministry of Education, Culture, Sports, Science and Technology (MEXT), Arctic Challenge for Sustainability project ArCS II (program grant number JPMXD1420318865).

# References
