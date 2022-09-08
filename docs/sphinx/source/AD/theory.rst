.. _theory:

Some FAQs
*********

Why use the Finite Differences when the adjoint is much faster?
===============================================================

The finite differences is supposed to be a validation test for the gradient generated using the adjoint mode. For fine grids, it only makes sense to use FD to validate the values of the adjoint at a few points since FD is expensive, requiring 2 non-linear forward model runs per evaluation of one component of the gradient. 
