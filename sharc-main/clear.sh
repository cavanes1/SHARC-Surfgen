#/bin/bash

# $PRIMARYDIR and SCRADIR are set
  PRIMARYDIR=$(pwd)
  TRAJ_DIR=$PRIMARYDIR/../trajdata-field
  rm -r $TRAJ_DIR/tmp*
  rm -r $PRIMARYDIR/outputdata/*
