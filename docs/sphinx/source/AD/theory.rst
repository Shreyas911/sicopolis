.. _theory:

Some FAQs
*********

Does the Finite Differences use AD at all?
==========================================

The Finite Differences does not use AD at all, it simply relies on multiple non-linear forward model runs to calculate the gradient. Thus Tapenade is not even required to differentiate the code. To activate the finite differences capabilities, the ``ALLOW_GRDCHK`` flag must be set, either within the compilation command itself or in the header file.



Why use the Finite Differences when the adjoint is much faster?
===============================================================

The finite differences is supposed to be a validation test for the gradient generated using the adjoint mode. For fine grids, it only makes sense to use FD to validate the values of the adjoint at a few points since FD is expensive, requiring 2 non-linear forward model runs per evaluation of one component of the gradient. 
