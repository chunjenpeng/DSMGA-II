#!/bin/bash
if [ "$#" != "2" ]; then
    echo "./run-nks1 FOLDER totalNP"
    exit 1
fi
FOLDER=$1
totalNP=$2
JOBNAME=NKs1
DATAPATH=/home/$USER/DSMGA-II/data
DIR=$DATAPATH/$FOLDER/$JOBNAME
mkdir -p $DIR

for ELL in 50 100 200 400
do
    for NUM in $(seq 0 99)
    do
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
        while [ $JOBNUM -ge $totalNP ]
        do
            sleep 1
            JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
        done
        ./sweep $ELL 10 4 1 $NUM > $DIR/$ELL-$NUM &
        echo "Submitting NKs1 $ELL-$NUM"
        sleep 1
    done
done

