#!/bin/bash

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# make_smb_corr.sh   Date: 2024-07-11
#
# -----------
# Description
# -----------
#
#   Extract surface mass balance (SMB) correction for a given time
#   from a SICOPOLIS output file.
#
# -------
# Example
# -------
#
#   ./make_smb_corr.sh \
#      grl05_spinup 0001 \
#      /work/sicopolis/grl05_spinup /home/sicopolis/sico_in/grl
#
#     -> create SMB correction (grl05_spinup0001_smb_corr.nc)
#           from output file grl05_spinup0001.nc
#              in directory /work/sicopolis/grl05_spinup,
#                 write it to directory /home/sicopolis/sico_in/grl
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

echo "Making smb_corr file..."

#-------- Settings --------

runname=$1
ergnum=$2
path1=$3
path2=$4

filename1="${path1}/${runname}${ergnum}.nc"
filename2="${path2}/${runname}${ergnum}_smb_corr.nc"

echo "   filename1 = ${filename1}"
echo "   filename2 = ${filename2}"

#-------- Extracting SMB correction from SICOPOLIS output file --------

ncks -O -F -v smb_corr ${filename1} ${filename2}
ncrename -v smb_corr,DSMB ${filename2}

#-------- End of script --------

echo "Done."

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
