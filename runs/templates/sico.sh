#!/bin/bash
LANG=C

################################################################################
#
#  s i c o . s h
#
#  bash script for
#  compilation, linking and execution of the program SICOPOLIS.
#
#  Authors: Malte Thoma, Thomas Goelles, Ralf Greve, Fuyuki Saito
#
#  Date: 2020-07-02
#
#    Execute script 
#       ./sico.sh -m <run_name> [further options...]
#    where <run_name> is the name of the simulation.
#
################################################################################

function error()
{
   echo -e "\n******************************************************************************">&2
   echo -e   "** ERROR:                                                                   **">&2
   echo -e   "$1" | awk '{printf("** %-73s**\n",$0)}'>&2
   echo -e "******************************************************************************\n">&2
   exit 1
}

################################################################################

function info()
{
   echo -e "$1" >&2
}

################################################################################

function usage()
{
   echo -e "\n  Usage: `basename $0` -m <run_name> [further options...]\n"\
   "     use -m? to get a list of available simulations\n"\
   "     [-m <run_name>] => name of the simulation\n"\
   "     [-i <dir>] => individual input directory, default is sico_in\n"\
   "     [-d <dir>] => individual output directory, default is sico_out/<run_name>\n"\
   "     [-a <ANFDATPATH>] => the ANFDATPATH directory, only needed if ANF_DAT=3,\n"\
   "                                         or if ANF_DAT=1 and TEMP_INIT=5\n"\
   "     [-t <TARGET_TOPO_DAT_PATH>] => the TARGET_TOPO_DAT_PATH directory,\n"\
   "                                    only needed if THK_EVOL=2 or 3,\n"\
   "                                    or if ACCSURFACE=7 and ABLSURFACE=7\n"\
   "     [-f] => force overwriting the output directory\n"\
   "     [-n] => skip make clean\n"\
   "     [-b] => skip execution, build only\n"\
   "     [-c <FILE>] => configuration file FILE instead of sico_configs.sh\n"\
   "     [-u] => (redundant option)\n"\
   "     [-z] => (redundant option)\n"\

      if [ "$1" ]; then error "$1"; fi
}

################################################################################

function check_args()
{
   while getopts a:bc:d:fhi:m:nt:uz? OPT ; do
     case $OPT in
       a) ANFDATPATH=$OPTARG ;;
       b) BUILD_ONLY="TRUE";;
       c) CONFIG=$OPTARG ;;
       d) local OUTDIR=$OPTARG ;;
       f) local FORCE="TRUE";;
       i) local INDIRIN=$OPTARG ;;
       m) RUN=$OPTARG ;;
       n) SKIP_MAKECLEAN="TRUE";;
       t) TARGET_TOPO_DAT_PATH=$OPTARG ;;
       u) local REDUNDANT_OPTION_U="TRUE";;
       z) local REDUNDANT_OPTION_Z="TRUE";;
     h|?) usage; exit 1;;
     
     esac            
   done

   PROGNAME="sicopolis"

   # Check if RUN is set correctly
   if [ ! "$RUN" ]; then error "No simulation set. Try option -h."; fi
   if [ $RUN == "?" ]; then
      info "--------------------------" 
      info "Available simulations are:"
      info "--------------------------"
      ls headers/ | sed 's/sico_specs_\(.*\)\.h/\1/'
      exit 1
   fi
   
   HEADER=`pwd`"/headers/sico_specs_$RUN.h"
   if [ ! -e $HEADER ]; then 
      error "Simulation header $HEADER does not exist."
   fi

   RUN_SPECS_HEADER="sico_specs_$RUN.h"

   # Output directory, absolute paths
   if [ ! "$OUTDIR" ]; then
      RESDIR=${PWD}"/../sico_out/$RUN"
   else
      lastch=`echo $OUTDIR | sed -e 's/\(^.*\)\(.$\)/\2/'`
      if [ ${lastch} == "/" ]; then OUTDIR=`echo $OUTDIR | sed '$s/.$//'`; fi
      if [ ! -e $OUTDIR ]; then error "$OUTDIR does not exist."; fi
      RESDIR=${OUTDIR}"/$RUN"
   fi
   
   # Input directory, absolute paths
      if [ ! "$INDIRIN" ]; then
      INDIR=${PWD}"/../sico_in"
   else
      lastch=`echo $INDIRIN | sed -e 's/\(^.*\)\(.$\)/\2/'`
      if [ ${lastch} == "/" ]; then INDIRIN=`echo $INDIRIN  | sed '$s/.$//'`; fi
      if [ ! -e $INDIRIN  ]; then error "$INDIRIN does not exist."; fi
      INDIR=${INDIRIN}
   fi
   
   # Checking for existing output
   if [ "$FORCE" ]; then $RM -rf $RESDIR 2> /dev/null ; fi
   if [ -e $RESDIR ]; then error "$RESDIR exists. Use -f to overwrite."; fi
   
   # Handling the ANFDATPATH
   # Reading variables from header
   ANF_DAT=$(sed -n 's%#define ANF_DAT % %p'  $HEADER)
   TEMP_INIT=$(sed -n 's%#define TEMP_INIT % %p'  $HEADER)
   ANFDATNAME=$(sed -n 's%#define ANFDATNAME % %p'  $HEADER | 
   sed -e "s%'%%g" | sed -e 's,^ *,,')
   SED_FOR_A="FALSE"
   if [[ $ANF_DAT -eq 3 || ($ANF_DAT -eq 1 && $TEMP_INIT -eq 5) ]]; then
      if [ ! "$ANFDATPATH" ]; then 
         error "ANFDATPATH not set. Use -a <ANFDATPATH>."; fi
         lastch=`echo $ANFDATPATH | sed -e 's/\(^.*\)\(.$\)/\2/'`
      # Handling of / at the end   
      if [ ${lastch} == "/" ]; then 
         ANFDATPATH=`echo $ANFDATPATH  | sed '$s/.$//'`
      fi
      # Check if file exists
      if [ ! -e ${ANFDATPATH}/${ANFDATNAME} ]; then 
         error "$ANFDATPATH/$ANFDATNAME does not exist."
      else
         SED_FOR_A="TRUE"       
      fi          
   fi
   
   # Handling of TARGET_TOPO_DAT_PATH
   THK_EVOL=$(sed -n 's%#define THK_EVOL % %p'  $HEADER)
   ACCSURFACE=$(sed -n 's%#define ACCSURFACE % %p'  $HEADER)
   ABLSURFACE=$(sed -n 's%#define ABLSURFACE % %p'  $HEADER)
   TARGET_TOPO_DAT_NAME=$(sed -n 's%#define TARGET_TOPO_DAT_NAME % %p' $HEADER | 
   sed -e "s%'%%g" | sed -e 's,^ *,,')
   SED_FOR_T="FALSE"
   if [[ ($THK_EVOL -eq 2 || $THK_EVOL -eq 3) || (ACCSURFACE -eq 7 && ABLSURFACE -eq 7) ]]; then
      if [ ! "$TARGET_TOPO_DAT_PATH" ]; then 
         error "TARGET_TOPO_DAT_PATH not set. Use -t <TARGET_TOPO_DAT_PATH>."; fi
      lastch=`echo $TARGET_TOPO_DAT_PATH | sed -e 's/\(^.*\)\(.$\)/\2/'`
      # Handling of / at the end   
      if [ ${lastch} == "/" ]; then 
         TARGET_TOPO_DAT_PATH=`echo $TARGET_TOPO_DAT_PATH  | sed '$s/.$//'`
      fi
      # Check if file exists
      if [ ! -e ${TARGET_TOPO_DAT_PATH}/${TARGET_TOPO_DAT_NAME} ]; then 
         error "$TARGET_TOPO_DAT_PATH/$TARGET_TOPO_DAT_NAME does not exist."
      else
         SED_FOR_T="TRUE"   
      fi         
   fi  

   # Configuration file
   if [ ! "$CONFIG" ]; then
      CONFIG="./sico_configs.sh"
   else
      if [ ! -e $CONFIG ]; then error "$CONFIG does not exist."; fi
      CONFIG=`dirname $CONFIG`/`basename $CONFIG`
   fi

   # Redundant options
   if [ $REDUNDANT_OPTION_U ]; then info "\nOption -u is redundant."; fi
   if [ $REDUNDANT_OPTION_Z ]; then info "\nOption -z is redundant."; fi
}

################################################################################

function compile()
{
   source ./sico_environment.sh

   info "\nConfiguration file ${CONFIG}."
   source $CONFIG

   cd ../src/

   EXE_FILE='sico_'${RUN}'.x'

   $RM -f ${EXE_FILE}
   $RM -f *.mod
   
   $CP $HEADER $RUN_SPECS_HEADER
   # Reading the header and changing the values
   sed 's%INPATH .*%INPATH '\'''${INDIR}\''%' $RUN_SPECS_HEADER |
   sed 's%OUTPATH .*%OUTPATH '\'''${RESDIR}\''%' > $$
   $MV $$ $RUN_SPECS_HEADER
   
   if [ $SED_FOR_A == "TRUE" ]; then
      sed 's%ANFDATPATH .*%ANFDATPATH\  '\'''${ANFDATPATH}\''%' $RUN_SPECS_HEADER > $$
      $MV $$ $RUN_SPECS_HEADER
   fi   
   
   if [ $SED_FOR_T == "TRUE" ]; then
      sed 's%TARGET_TOPO_DAT_PATH .*%TARGET_TOPO_DAT_PATH\  '\'''${TARGET_TOPO_DAT_PATH}\''%' $RUN_SPECS_HEADER > $$
      $MV $$ $RUN_SPECS_HEADER
   fi

   FCFLAGS="${FCFLAGS} -DRUN_SPECS_HEADER=\"${RUN_SPECS_HEADER}\" -o ${EXE_FILE}"
   ${FC} ./${PROGNAME}.F90 ${FCFLAGS}

   mkdir $RESDIR
   
   $MV $RUN_SPECS_HEADER $RESDIR
   
   ### # Writing a log file with some information about Subversion
   ### SVNINFOFILE=$RESDIR/svninfo.log
   ### echo "Subversion infos:" > $SVNINFOFILE
   ### echo "-----------------" >> $SVNINFOFILE
   ### echo -n "Revision no. of src directory: " >> $SVNINFOFILE
   ### svnversion >> $SVNINFOFILE
   ### svn info ${PROGNAME}.F90 | grep "Repository Root" >> $SVNINFOFILE
   
   # Writing a log file with some information about the host
   HOSTINFOFILE=$RESDIR/hostinfo.log
   echo "Host infos:" > $HOSTINFOFILE
   echo "-----------" >> $HOSTINFOFILE
   echo -n "Host name: " >> $HOSTINFOFILE
   hostname >> $HOSTINFOFILE
   echo -n "OS: " >> $HOSTINFOFILE
   uname >> $HOSTINFOFILE
   echo -n "User: " >> $HOSTINFOFILE
   id >> $HOSTINFOFILE

   cd $OLDPWD
}

################################################################################

function run()
{
   if [ $BUILD_ONLY ]; then
      info "Skip execution."
      return 0
   fi

   cd ../src/

   # Needed for openmp on LINUX 
   UNAME=`uname`
   if [ "$UNAME" == "Linux" ]; then
      echo "Setting ulimit to unlimited."
      ulimit -s unlimited
      export STACKSIZE=8192
   fi
   
   local OUT=../runs/out_$RUN.dat
   info "Starting ./$EXE_FILE"
   info "         (log-output in out_$RUN.dat) ..."

   (time $NICE -n 19 ./${EXE_FILE}) >$OUT
   # (time $NICE -n 19 ./${EXE_FILE}) >$OUT 2>&1
   
   $MV $OUT $RESDIR
  
   if [ ! $SKIP_MAKECLEAN ]; then
      $RM -f ./${EXE_FILE}
      $RM -f *.mod
   fi

   info "\n... finished"
   info   "    (output in $RESDIR).\n"

   cd $OLDPWD
}

################################################################################

RM=/bin/rm
CP=/bin/cp
MV=/bin/mv

if [ -x /usr/bin/nice ]; then
   NICE=/usr/bin/nice
elif [ -x /bin/nice ]; then
   NICE=/bin/nice
else
   NICE=nice
fi

check_args $*
compile 
run

############################ End of sico.sh ####################################
