#!/bin/bash
if [ "$#" != "2" ]; then
    echo "./run-other FOLDER totalNP"
    exit 1
fi
FOLDER=$1
totalNP=$2
DATAPATH=/home/$USER/DSMGA-II/data

JOBNAME=ONEMAX
DIR=$DATAPATH/$FOLDER/$JOBNAME
NUM="0"
mkdir -p $DIR
for ELL in 100 200 400
do
    JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    while [ $JOBNUM -ge $totalNP ]
    do
        sleep 1
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    done
    ./sweep $ELL 10 0 $NUM > $DIR/$ELL-$NUM &
    echo "Submitting $JOBNAME-$ELL"
    sleep 1
done

JOBNAME=MK
DIR=$DATAPATH/$FOLDER/$JOBNAME
NUM="0"
mkdir -p $DIR
for ELL in 100 200 400
do
    JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    while [ $JOBNUM -ge $totalNP ]
    do
        sleep 1
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    done
    ./sweep $ELL 10 1 $NUM > $DIR/$ELL-$NUM &
    echo "Submitting $JOBNAME-$ELL"
    sleep 1
done

JOBNAME=FTRAP
DIR=$DATAPATH/$FOLDER/$JOBNAME
NUM="0"
mkdir -p $DIR
for ELL in 120 240 480
do
    JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    while [ $JOBNUM -ge $totalNP ]
    do
        sleep 1
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    done
    ./sweep $ELL 10 2 $NUM > $DIR/$ELL-$NUM &
    echo "Submitting $JOBNAME-$ELL"
    sleep 1
done

JOBNAME=CYC
DIR=$DATAPATH/$FOLDER/$JOBNAME
NUM="0"
mkdir -p $DIR
for ELL in 100 200 400
do
    JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    while [ $JOBNUM -ge $totalNP ]
    do
        sleep 1
        JOBNUM=$(ps -ef | grep $USER | grep sweep | wc -l)
    done
    ./sweep $ELL 10 3 $NUM > $DIR/$ELL-$NUM &
    echo "Submitting $JOBNAME-$ELL"
    sleep 1
done
