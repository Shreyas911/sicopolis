#!/bin/bash

function compare_gradients_with_fd(){

	if [ "$#" -ne 3 ]; then
		echo "compare_outputs Usage: Need two data files and 1 string as arguments. First file - FD, Second file - adjoint/forward. First string - forward/adjoint."
		exit 1
	fi

	FD_FILE="$1"
	AD_FILE="$2"
	CTRL_FLOW="$3"

	if [ "$(cat ${FD_FILE} | wc -l)" -ne "$(cat ${AD_FILE} | wc -l)" ]; then
		echo "Different number of lines in both files."
		exit 1
	fi

	DIFF=$(diff -w <(head ${FD_FILE}) <(head ${AD_FILE}))

	if [ "${DIFF:-0}" == 0 ]; then
		echo "Adjoint and FD match exactly, highly unlikely."
		exit 1
	else
		#echo $(diff -w -y <(head ${FD_FILE}) <(head ${AD_FILE})) 
		readarray FD_VALS < $FD_FILE
		readarray AD_VALS < $AD_FILE	
	fi		
	
	NUM_POINTS=$(cat ${AD_FILE} | wc -l)
	for (( i = 0; i <= $NUM_POINTS-1; i++)) do
		GRDCHK_VAL=$(printf "${FD_VALS[i]}\n")
		AD_VAL=$(printf "${AD_VALS[i]}\n")
		RATIO=$(bc -l <<< '(100*('$GRDCHK_VAL'-'$AD_VAL')/'$GRDCHK_VAL')')
		echo ${RATIO}
		if [ $(echo "scale = 10; $RATIO < -10.0 " | bc -l) = 1 ] || [ $(echo "scale = 10; $RATIO > 10.0 " | bc -l) = 1 ]; then
			echo "Relative error greater than 5%"
			exit 1
		fi
	done
}
compare_gradients_with_fd $1 $2 $3
