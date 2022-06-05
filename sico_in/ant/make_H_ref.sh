#!/bin/bash

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# make_H_ref.sh   Date: 2022-06-05
#
# Description:
#
#   Generate reference ice thickness for a given time from a
#   SICOPOLIS output file.
#
#   Example:
#
#      ./make_H_ref.sh \
#         -d /work/sicopolis/ant08_hist \
#         -m ant08_hist -n 0001 -y 2015
#
#      -> create reference ice thickness (ant08_hist0001_H_ref_2015.nc)
#            from output file ant08_hist0001.nc
#               in directory /work/sicopolis/ant08_hist
#                  for the year 2015.
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# set -x

#-------- Flags --------

while getopts d:m:n:ty:z: flag
do
   case ${flag} in
      d) outdir=${OPTARG};;
      m) runname=${OPTARG};;
      n) ergnum=${OPTARG};;
      t) two_d_file="TRUE";;
      y) year=${OPTARG};;
   esac            
done

if [ $two_d_file ]; then
   two_d="_2d_"
else
   two_d=""
fi

#-------- Settings --------

dir1="${outdir}"
dir2="."

filename1="${runname}${two_d}${ergnum}.nc"
filename2="${runname}${two_d}${ergnum}_H_ref_${year}.nc"
filename3="${dir1}/${filename1}"
filename4="${dir2}/${filename2}"

#-------- Extracting ice thickness from SICOPOLIS output file --------

echo "Creating ${filename2}..." ;

# ncks -O -F -v crs,x,y,lon,lat,H ${filename3} ${filename4}
ncks -O -F -v mapping,x,y,lon,lat,H ${filename3} ${filename4}

#-------- End of script --------

echo "Done." ;

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
