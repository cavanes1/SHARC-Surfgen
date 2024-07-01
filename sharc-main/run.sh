#/bin/bash

# $PRIMARYDIR and SCRADIR are set
ST=${1:-1}
EN=${2:-10}
for ((i = ST; i <= EN; i++))
do 
  PRIMARYDIR=$(pwd)
  #DATADIR=$PRIMARYDIR/outputdata
  COPY_DIR=$PRIMARYDIR/../trajdata/traj$i
  cd $COPY_DIR
  $SHARC/sharc.x input
  err=$?
  if [ ! $err == 0 ];
    then
      echo "abnormal abort"
   #  cp $COPY_DIR/QM/* $PRIMARYDIR/QM/
  fi
#  $SHARC/data_extractor.x -e output.dat

#  cp $COPY_DIR/output.lis $DATADIR/output.lis.$i
#  cp $COPY_DIR/output_data/energy.out $DATADIR/energy.out.$i
#  rm -r $COPY_DIR
#  exit $err
  cd $PRIMARYDIR
  echo "finish traj", $i
done
