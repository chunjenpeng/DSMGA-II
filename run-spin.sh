#!/bin/bash
if [ "$#" != "2" ]; then
    echo "./run-spin FOLDER totalNP"
    exit 1
fi
FOLDER=$1
totalNP=$2
JOBNAME=SPIN
DATAPATH=/home/$USER/DSMGA-II/data
DIR=$DATAPATH/$FOLDER/$JOBNAME
mkdir -p $DIR

for ELL in 36 100 196 400 
do
    for NUM in $(seq 1 100)
    do
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
        while [ $JOBNUM -ge $totalNP ]
        do
            sleep 1
            JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
        done
        ./sweep $ELL 10 5 $NUM > $DIR/$ELL-$NUM &
        echo "Submitting SPIN $ELL-$NUM"
        sleep 1
    done
done

