#!/bin/bash

echo "==> Making smb_corr file ..." ;

runname=$1
ergnum=$2
path1=$3
path2=$4

filename1="${path1}/${runname}${ergnum}.nc"
filename2="${path2}/${runname}${ergnum}_smb_corr.nc"

echo "==> filename1 = ${filename1}"
echo "==> filename2 = ${filename2}"

ncks -O -F -v smb_corr ${filename1} ${filename2}
ncrename -v smb_corr,DSMB ${filename2}

echo "==> Done." ;
#
