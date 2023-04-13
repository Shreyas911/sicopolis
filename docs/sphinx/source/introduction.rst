.. _introduction:

Introduction
************

SICOPOLIS (SImulation COde for POLythermal Ice Sheets) is a 3-d dynamic/thermodynamic model that simulates the evolution of large ice sheets and ice caps. It was originally created by Greve (1997a,b) in a version for the Greenland ice sheet. Since then, SICOPOLIS has been developed continuously and applied to problems of past, present and future glaciation of Greenland, Antarctica, the entire northern hemisphere, the polar ice caps of the planet Mars and others.

The model is based on the shallow ice approximation for grounded ice, the shallow shelf approximation for floating ice (e.g., Greve and Blatter 2009) and, optionally, hybrid shallow-ice–shelfy-stream dynamics for ice streams (Bernales et al. 2017a,b). It is coded in Fortran and uses finite difference discretisation on a staggered (Arakawa C) grid, the velocity components being taken between grid points. A variety of different thermodynamics solvers are available, namely the polythermal two-layer method, two versions of the one-layer enthalpy method, the cold-ice method and the isothermal method (Greve and Blatter 2016).

The coding is based on a consequent low-tech philosophy. All structures are kept as simple as possible, and advanced coding techniques are only employed where it is deemed appropriate. The use of external libraries is kept at an absolute minimum, which makes the installation very easy and fast.

Legal notes
===========

SICOPOLIS is free and open-source software. It can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at the user's option) any later version.

SICOPOLIS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
