#!/bin/bash
#
# SetCVSRoot.sh: This script will replace all CVS/Root files with
#   the value of the first argument to this script.
# ie: share/Scripts/SetCVSRoot.sh darrens@herot:/FRAMEWORK

echo Setting CVS/Root to $1 in $PWD

FILES=`find . -type f -name Root`

for file in $FILES
do
  echo Fixing $file
  echo $1 > $file
done

exit 0
