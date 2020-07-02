#!/bin/bash
LANG=C

################################################################################
#
#  r e v _ i d . s h
#
#  bash script for determining a unique revision identifier from Git.
#
#  Author: Ralf Greve
#
#  Date: 2020-07-03
#
################################################################################

git status >/dev/null 2>&1

if [[ `echo $?` -eq 0 ]]; then

   REPO_REVISION_=`git rev-list HEAD --count`
   BUILD_BRANCH=`git rev-parse --abbrev-ref HEAD`
   BUILD_REV_ID=`git rev-parse HEAD`
   BUILD_REV_ID_SHORT=`git describe --long --tags --dirty --always`
   if [ ${BUILD_BRANCH} == "master" ]; then
      REPO_REVISION=${REPO_REVISION_}_g${BUILD_REV_ID_SHORT}
   else
      REPO_REVISION=${BUILD_BRANCH}_${REPO_REVISION_}_r${BUILD_REV_ID_SHORT}
   fi
   ### Source: https://stackoverflow.com/a/29922075

else

   REPO_REVISION="(no Git repository)"
   BUILD_BRANCH=
   BUILD_REV_ID=

fi

export REPO_REVISION
export BUILD_BRANCH
export BUILD_REV_ID

echo "Git revision identifier: ${REPO_REVISION}"

########################### End of rev_id.sh ###################################
