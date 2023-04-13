.. _introduction:

Introduction
************

SICOPOLIS (SImulation COde for POLythermal Ice Sheets) is a 3-d dynamic/thermodynamic model that simulates the evolution of large ice sheets and ice caps. It was originally created by Greve :cite:`greve_1997a, greve_1997b` in a version for the Greenland ice sheet. Since then, SICOPOLIS has been developed continuously and applied to problems of past, present and future glaciation of Greenland, Antarctica, the entire northern hemisphere, the polar ice caps of the planet Mars and others.

The model is based on the shallow ice approximation for grounded ice, the shallow shelf approximation for floating ice (e.g., Greve and Blatter :cite:`greve_blatter_2009`) and, optionally, hybrid shallow-ice–shelfy-stream dynamics for ice streams (Bernales et al. :cite:`bernales_etal_2017a, bernales_etal_2017b`). It is coded in Fortran and uses finite difference discretisation on a staggered Arakawa C grid (Arakawa and Lamb :cite:`arakawa_lamb_1977`), the velocity components being taken between grid points. A variety of different thermodynamics solvers are available, namely the polythermal two-layer method, two versions of the one-layer enthalpy method, the cold-ice method and the isothermal method (Greve and Blatter :cite:`greve_blatter_2016`).

The coding is based on a consequent low-tech philosophy. All structures are kept as simple as possible, and advanced coding techniques are only employed where it is deemed appropriate. The use of external libraries is kept at an absolute minimum, which makes the installation very easy and fast.

Legal notes
===========

| Copyright 2009–2023 SICOPOLIS Authors
| (Ralf Greve, Jorge Bernales, Sebastian Beyer, Heinz Blatter, Reinhard Calov, Thorben Dunse, Eduardo Flandez, Shreyas Sunil Gaikwad, Ben Galton-Fenzi, Thomas Gölles, Björn Grieger, Philipp Hancke, Laurent Hascoet, Patrick Heimbach, Nina Kirchner, Thomas Kleiner, Sascha Knell, Anne Le Brocq, Liz Curry-Logan, Sri Hari Krishna Narayanan, Alex Robinson, Fuyuki Saito, Tatsuru Sato, Marius Schäfer, Matthias Scheiter, Oliver J. Stenzel, Malte Thoma, Roland Warner).

SICOPOLIS is free and open-source software. It can be redistributed and/or modified under the terms of the `GNU General Public License <https://www.gnu.org/licenses/>`__ as published by the Free Software Foundation, either version 3 of the License, or (at the user's option) any later version.

SICOPOLIS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the `GNU General Public License <https://www.gnu.org/licenses/>`__ for more details.
