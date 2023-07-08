.. _ad_intro_features:

Introduction and features
=========================

Previously, OpenAD have been used to get the adjoint of the SICOPOLIS code (Logan et al. :cite:`logan_etal_2020`). The current implementation with Tapenade, SICOPOLIS-AD v2 (Gaikwad et al. :cite:`gaikwad_etal_2023`), has the following advantages over the previous implementation\:

1. It is up-to-date with the latest SICOPOLIS code

2. The AD tool Tapenade is open-source and actively maintained

3. A new tangent linear code generation capability is introduced (Forward Mode)

   * This is useful for adjoint validation, and from a physical perspective, also useful for Bayesian UQ and inverse modeling.
  
4. We are now able to deal with inputs in the NetCDF format

5. We have now correctly incorporated the external LIS solver, its tangent linear code, and its adjoint which improve the simulation of Antarctic ice shelves and Greenland outlet glaciers. (subject to more testing)

6. We leverage continuous integration and the pytest framework in order to track changes in the trunk that "break" the AD-based code generation - precludes the need for constant monitoring.

7. We "show" the entire code to Tapenade, including the initialization subroutines, thus avoiding cumbersome maintenance of subroutines OpenAD used to initialize for adjoint runs.

8. We have provided convenient Python scripts to make I/O with the differentiated variables easier.

In addition we also have the following previously available features\: 

1. The adjoint mode is available, like before along with the capability to do Finite Differences validation of the gradient computed using the adjoint mode.

The code has the following capabilities or possible applications\: 

1. Paleoclimatic inversions using the adjoint generated gradients as part of a model calibration exercise.

2. Uncertainty Quantification of calibrated parameters by leveraging the new tangent linear mode and adjoint mode.

3. Sensitivity analysis of various quantities of interest (QoI) to state parameters.

4. Optimal Experimental Design - where should sensors be placed such that the newly collected data optimally informs our uncertain parameters.
