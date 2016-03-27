#!/bin/bash
if [ "$#" == "0" ]; then
     echo "./average DIR ell"
     exit 1
fi

DIR=$1
ELL=$2

grep NFE $DIR/$ELL-* > tempFile

count=0
total1=0
total2=0

awk '{ total1 += $2; count++} END { print "average NFE =", total1/count}' tempFile
#
#grep ls $DIR/$ELL-* > tempFile
#
#count=0
#total1=0
#total2=0
#
#awk '{ total1 += $2; total2 += sqrt($3/9); count++ } END { print 
#total1/count; print total2/count }' tempFile

wc -l tempFile
