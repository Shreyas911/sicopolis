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
#  Date: 2023-03-03
#
################################################################################

git status >/dev/null 2>&1

if [[ `echo $?` -eq 0 ]]; then

   # REPO_REVISION_=`git rev-list HEAD --count`
   BUILD_BRANCH=`git rev-parse --abbrev-ref HEAD`
   BUILD_REV_ID=`git rev-parse HEAD`
   BUILD_REV_ID_SHORT=`git rev-parse --short=8 HEAD`

   REPO_REVISION="${BUILD_BRANCH}_${BUILD_REV_ID_SHORT}"

   if [[ -n $(git status --untracked-files=no --porcelain) ]]; then
      REPO_REVISION="${REPO_REVISION}_dirty"
   fi

   ### Source: https://stackoverflow.com/a/29922075 (modified)

else

   REPO_REVISION="no_Git_repository"

fi

export REPO_REVISION

echo "Git revision identifier: ${REPO_REVISION}"

########################### End of rev_id.sh ###################################
