#!/usr/bin/bash

if [ $1 == "-h" -o $1 == "--help" -o $# == 0 ]
then printf "\n
Usage:
    sh $0 column_number times 'start'
    sh $0 column_number times 'pos'
    
if choose 'start':
    sh $0 5 4 start
    will give: 1,2,3,4,5,10,15,20

if choose 'pos':
    sh $0 5 4 pos
    will give: 5,10,15,20

"
     exit 0
fi


if [ $3 == "start" ]
then
    echo ""| awk -v col=$1 -v times=$2 'BEGIN{OFS=",";ORS=","}{for(i=1;i<col;i++){print i} for(i=1;i<times;i++){print col*i} if(i==times){printf col*i}}'
fi


if [ $3 == "pos" ]
then
    echo ""| awk -v col=$1 -v times=$2 'BEGIN{OFS=",";ORS=","}{for(i=1;i<times;i++){print col*i} if(i==times){printf col*i}}'
fi


