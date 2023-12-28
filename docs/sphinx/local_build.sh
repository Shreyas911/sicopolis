#!/bin/bash
# (Selection of shell)

#--------------------------------------------------------------------------
# Building the SICOPOLIS manual with Sphinx locally.
#--------------------------------------------------------------------------

# Requirements:
# -------------
# Ensure you have pip (package installer for Python) installed
# (https://pypi.org/project/pip/).

# cd docs/sphinx
# pip install -r requirements.txt

# Actual build:
# -------------

make clean
make html
# make latexpdf PAPER=a4

# Open docs/sphinx/build/html/index.html to see the locally built manual.

#--------------------------------------------------------------------------
#
