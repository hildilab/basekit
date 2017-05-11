#!/bin/bash

fname=$(basename "$1")

cp -rf /home/student/Software/Dowser++/files $2
fname=$(python /home/student/Documents/basekit/basekit/parameterizationDowserPlus2.py $fname $2)

cd $2/files
cp $2/$fname $2/files
./main.sh $fname $3
cp PredictedInternal.pdb ..
cp PredictedInternal.pdb ../Dowser++_Water