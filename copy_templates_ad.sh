#!/bin/bash
# (Selection of shell)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   Making local, not version-controlled copies of the scripts and the
#   run-specification header files required for running SICOPOLIS-AD.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CP=/bin/cp

RUN_DIR=./runs
HEADER_DIR=${RUN_DIR}/headers

SRC_DIR=./src
AD_DIR=${SRC_DIR}/subroutines/tapenade

echo "Copying scripts"
echo "        from ${AD_DIR} to ${SRC_DIR} ..."
echo " "
$CP -f ${AD_DIR}/regression_test.sh ${SRC_DIR}
$CP -f ${AD_DIR}/preprocessor.py ${SRC_DIR}

echo "Copying run-specification headers"
echo "        from ${AD_DIR} to ${HEADER_DIR} ..."
echo " "
$CP -f ${AD_DIR}/sico_specs_*.h ${HEADER_DIR}

echo "... done!"
echo " "

#--------
