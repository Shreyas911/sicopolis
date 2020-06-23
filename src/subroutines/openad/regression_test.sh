#!/bin/bash
#__________________________________________________________________________________
# This script is to be run to ensure that the adjoint code of SICOPOLIS is properly
# functioning. That entails comparing two types of outputs: (1) costs (e.g., the 
# volume of an ice sheet, or model misfit to some data type, like ages); and (2) 
# gradients (i.e., change in cost due to some perturbation to a control parameter 
# -- termed "finite differences" or "gradient check" in the below).
# Gradients may be generated via two different methods in SICOPOLIS: (a) finite-difference
# approximations, or the gradient check below that employs only the original, forward
# SICOPOLIS code; and (b) by the adjoint code, which is compiled and executed herein. 

# For a detailed explanation of the adjoint code structure (e.g., file placement, 
# choice of control variables and cost function) please refer to the Sicopolis-AD
# Quickstart and User manual (https://www.overleaf.com/read/xtybnbhxmbcc).

# Last modified: Liz Logan, 4/5/2018
#__________________________________________________________________________________

#__________________________________________________________________________________
function page_break()
# A cosmetic function that indicates to the user when significant portions of this
# regression_test have been completed (note: not necessarily successfully). Simply 
# outputs to screen lines and statement that that portion has completed.
#__________________________________________________________________________________
{
HEADER_LOCAL="$1"

for line in {1..15}
do
   echo "-------------------------------------------------------------------------"
done
echo "                    "$HEADER_LOCAL" portion complete."   
}

#__________________________________________________________________________________
function compare_outputs()
# Compares regression_test results for the gradient check vs. adjoint values. This 
# function is sensitive to output file location and number of points you have chosen
# to compare. If in file src/subroutines/openad/opeand_m.F90 you modified the number
# of gradient check points from the default 10 to some number NUM, please change
# that below. If some error occurs here, it is 
# likely due to (1) - model did not output results; (2) - user has modified directory
# output location or structure.
#__________________________________________________________________________________
{
# If you changed the number of points to be compared between FD and AD modes, please
# indicate that number here:
NUM_POINTS=10
# Non-clever solution for bash counting from 0:
let NUM_POINTS--

# reassigning input arguments
# options (currently) are:
VALUES="$1"            # COSTS, GRADIENTS
HEADER_LOCAL="$2"      # this comes from the HEADER_FILE options in the main loop below 
CTRL_FLOW="$3"         # forward, adjoint
MODE="$4"              # (adjoint only, modes plain, tape, adjoint): -p, -t, -a
OUT_PATH_LOCAL="$5"    # output path 

# comparing the costs:
if [ $VALUES == "COSTS" ]; then

# Grab the gradient check cost files:
  FILE_GRDCHK=$OUT_PATH_LOCAL"/CostVals_"$HEADER_LOCAL"_GRDCHK.dat"
  COST_GRDCHK=$(head -n 1 $FILE_GRDCHK)  

# Grab the adjoint-modes cost files:
  FILE_OAD=$OUT_PATH_LOCAL"/AD_COST_"$HEADER_LOCAL"_OAD"$MODE".dat"
  COST_OAD=$(head -n 1 $FILE_OAD)

# Calculate the ratio 100*[Cost(gradient check) - Cost(adjoint)] / Cost(gradient check)
  RATIO=$(bc -l <<< '(100*('$COST_GRDCHK'-'$COST_OAD')/'$COST_GRDCHK')')

# Write this percent different to the regression_log:
  echo  $RATIO" percent different in "$MODE" mode."                              >> regression_log.txt

# comparing the gradients:
elif [ $VALUES == "GRADIENTS" ]; then

echo " ______________________________________________________________"           >> regression_log.txt
echo " "                                                                         >> regression_log.txt
echo " Gradient comparison for: "$HEADER_LOCAL                                   >> regression_log.txt
echo "          FINITE-DIFFERENCES:                    ADJOINT VALUES:"          >> regression_log.txt  
echo " "                                                                         >> regression_log.txt

# Grab the gradient values from grdchk and adjoint sensitivity values from adjoint: 
  FILE_GRDCHK=$OUT_PATH_LOCAL"/GradientVals_"$HEADER_LOCAL"_GRDCHK.dat"
  FILE_AD=$OUT_PATH_LOCAL"/AD_grdchk_"$HEADER_LOCAL".dat"

# Preparing to compare the two columns of values and check if they're different:
  DIFF=$(diff -w -y -W 80 <(head ${FILE_GRDCHK}) <(tail -n +2 ${FILE_AD}))
 
  readarray      GRDCHK_VALS  < $FILE_GRDCHK
  readarray -s 1 AD_VALS      < $FILE_AD

# Writing results of value comparison to regression_log:
  if [ "${DIFF:-0}" == 0 ]; then
    echo " FD and adjoint values match exactly (computationally impossible)."    >> regression_log.txt
    echo " Something went wrong: please contact liz.curry.logan@gmail.com with " >> regression_log.txt
    echo " header file and source code. "                                        >> regression_log.txt
  else
    echo "$DIFF"                                                                 >> regression_log.txt
    echo " "                                                                     >> regression_log.txt
    echo "--- percent different (FD - AD / FD) ---"                              >> regression_log.txt
    echo " "                                                                     >> regression_log.txt
    for (( i = 0; i <= $NUM_POINTS; i++)) do 
      GRDCHK_VAL=$(printf "${GRDCHK_VALS[i]\n}")
      AD_VAL=$(printf "${AD_VALS[i]}\n")
#      if [ "$GRDCHK_VAL" == "0.00000000000000000000" ]; then 
#      elif (( "$GRDCHK_VAL" == 0.0)) && (( "$AD_VAL" != 0.0)); then
#        RATIO=100.0
#      else 
        RATIO=$(bc -l <<< '(100*('$GRDCHK_VAL'-'$AD_VAL')/'$GRDCHK_VAL')')
#      fi
      echo $RATIO                                                                >> regression_log.txt
    done
  fi

fi

}

#__________________________________________________________________________________
#__________________________________________________________________________________
#__________________________________________________________________________________
#
#                Main regression testing routine 
#__________________________________________________________________________________
#__________________________________________________________________________________
#__________________________________________________________________________________

SRC_PATH=$PWD
OUT_PATH=$SRC_PATH/../sico_out/OAD

# Default header file options to be tested (if you wish to test a custom one, 
# modify the header string here): 
#declare -a HEADER_FILE=("ant20_shelves")
#declare -a HEADER_FILE=("simpleaf")
#declare -a HEADER_FILE=("new_ref_output_20km_vol")
declare -a HEADER_FILE=("new_ref_output_20km")

# Arguments passed to executable in adjoint mode 
# ("-p" == plain,
#  "-t" == tape,
#  "-a" == adjoint): 
declare -a MODES=("-p" "-t" "-a")

# Remove all old adjoint outputs: 
if [ -e          $OUT_PATH ]; then
   if [ "$(ls -A $OUT_PATH)" ]; then
      rm $OUT_PATH/*
      echo "SicoAD outputs removed."
   fi
else 
   mkdir $OUT_PATH
   echo "SicoAD created output dir OAD in: " 
   echo $OUT_PATH
fi

# move around the maths modules: 
cp $SRC_PATH/subroutines/openad/*_maths_* \
   $SRC_PATH/subroutines/general

# Remove old regression_log:
if [ -e regression_log.txt ]; then
   rm $SRC_PATH/regression_log.txt
   echo "Previous regression_test.sh results removed."
fi

#changing into src dir for all the execution below:
cd $SRC_PATH

# Starting the clock on the total execution time of this script:
SECONDS=0

#___________________________________________
# BELOW:
#
# Generating gradients via: 
# 1 - finite differences;
# 2 - adjoint code.
# And:
# 3 - comparing the associated costs and 
#     gradient values. 
#___________________________________________

# Looping over all the adjoint tests desired:
for test in "${HEADER_FILE[@]}"
do

#___________________________________________
#
# 1 - Gradient check:
#___________________________________________
   HEADER=$test"_GRDCHK" 
   echo "Regression test: "$HEADER

#  compile gradient check driver:
   make -f MakefileOpenAD clean
#  simulations using Library of Iterative Solvers (LIS), or not:
   if [ -z "${LISDIR}" ] ; then
      make -f MakefileOpenAD LISDIR=$LISDIR HEADER=$HEADER drivergrdchk
   else
      make -f MakefileOpenAD                HEADER=$HEADER drivergrdchk
   fi

#  call the executable something else:
   EXECUTABLE=$HEADER".exe"
   mv drivergrdchk $EXECUTABLE 

#  execute gradient check: 
   time ./$EXECUTABLE 

#  move its output to adjoint output directory:
   mv GradientVals_$HEADER.dat   $OUT_PATH
   mv CostVals_$HEADER.dat       $OUT_PATH

#  prettify 
   page_break $HEADER

#___________________________________________
#
# 2 - Adjoint:
#___________________________________________
   HEADER=$test"_OAD"

#  compiling adjoint code:
   make -f MakefileOpenAD clean
#  simulations using Library of Iterative Solvers (LIS), or not:
   if [ -z "${LISDIR}"]; then
      make -f MakefileOpenAD LISDIR=$LISDIR HEADER=$HEADER driver
   else
      make -f MakefileOpenAD                HEADER=$HEADER driver
   fi

#  call the executable something else:
   EXECUTABLE=$HEADER".exe"
   mv driver $EXECUTABLE 

##_____________________
##  loop over adjoint MODES: 
##       plain   (-p), 
##       tape    (-t), 
##       adjoint (-a)
##_____________________
   for mode in "${MODES[@]}"
   do

#  remove forward (gradient-checked) mode outputs:
      if [ -d  $OUT_PATH/../$HEADER ]; then
         rm -r $OUT_PATH/../$HEADER
      fi

#  file name for screen output:
      OUT_FNAME="out_"$HEADER$mode".dat"

#  execute adjoint code: 
      time ./$EXECUTABLE $mode > $OUT_FNAME 

#  move the output files (and re-write names):
      mv AD_COST    $OUT_PATH/"AD_COST_"$HEADER$mode".dat"
      mv $OUT_FNAME $OUT_PATH

#_____________________
#
# 3 - Compare and report: 
#   a - gradient check COSTS     vs. SicoAD 
#   b - gradient check GRADIENTS vs. SicoAD 
#_____________________
#   a - COSTS: 
      if [ $mode == "-p" ]; then 
          echo " ______________________________________________________________" >> regression_log.txt
          echo " "                                                               >> regression_log.txt
          echo " Cost comparison for:  "$HEADER >> regression_log.txt
          echo " ______________________________________________________________" >> regression_log.txt
      fi
      compare_outputs "COSTS" $test "adjoint" $mode $OUT_PATH

#   b - GRADIENTS: 
      if [ $mode == "-a" ]; then 
         mv AD_Vals_for_grdchk.dat    $OUT_PATH/"AD_grdchk_"$test".dat"
         compare_outputs "GRADIENTS" $test "adjoint" $mode $OUT_PATH
      fi
 
# prettify 
      page_break $HEADER$mode

   done
done

# Go back to hiding the maths modules in openad/ 
rm $SRC_PATH/subroutines/general/*_stub* 
rm $SRC_PATH/subroutines/general/*_grad*

# Stopping the clock on this total execution time:
ELAPSED="$(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

#cd back into cron dir to report results of regression
cd $CRON_PATH
