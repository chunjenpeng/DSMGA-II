#!/bin/bash
if [ "$#" != "2" ]; then
    echo "./run-sat.sh FOLDER totalNP"
    exit 1
fi
FOLDER=$1
totalNP=$2
JOBNAME=SAT
DATAPATH=/home/$USER/DSMGA-II/data
DIR=$DATAPATH/$FOLDER/$JOBNAME
mkdir -p $DIR

for ELL in 20 50 75 100 200 
do
    for NUM in $(seq 1 100)
    do
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
        while [ $JOBNUM -ge $totalNP ]
        do
            sleep 1
            JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
        done
        ./sweep $ELL 10 6 $NUM > $DIR/$ELL-$NUM &
        echo "Submitting SAT $ELL-$NUM"
        sleep 1
    done
done

