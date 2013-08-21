#!/bin/bash

DIST=`pwd`

SVMSRC=$DIST/src	# either manually set or use the loop below

PATCH=`which patch`		# location of patch command


if [ ! -x "$PATCH" ]; then
	echo "ERROR: could not find 'patch' command in your PATH environment."
fi

# if the SVMSRC directory isn't valid, ask where it is
while [ ! -f "$SVMSRC/svm_common.c" ];
do
  echo "Enter directory containing SVM v6 source: "
  read SVMSRC

  if [ ! -f "$SVMSRC/svm_common.c" ];
  then
	echo "Warning: did not find svm_common.c in $SVMSRC, try again."
  fi
done

if [ ! "$SVMSRC" == "$DIST/src" ]; then
   echo "Copying SVM-Lite files from $SVMSRC to $DIST/src"
   cp $SVMSRC/* $DIST/src
fi

echo "Copying distribution files from $DIST/dist to $DIST/src"
cp $DIST/dist/* $DIST/src

echo "Patching $SVMSRC"
cd $DIST/src

patch -p1 < $DIST/dist/patch.svm60

