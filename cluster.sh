#!/bin/bash

#cd to folder containing cluster.sh and then execute:
#./cluster.sh  full_path_to_folder  [start_sec]  [stop_sec]
#./cluster.sh  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/

FS=450450;
NW=15
K=29
PVAL=0.01
NFFT=(0.001 0.0005 0.00025)

for i in $(ls -1 $1*.ch* | sed s/.ch[0-9]// | uniq)
do
#  echo $i
  job_name=$(basename $i)
  for j in ${NFFT[@]}
  do
#    echo $j
    qsub -N "MTBP-$job_name-$j" -pe batch 8 -b y -j y -cwd -o "$job_name-$j.log" -V ./cluster2.sh "\"$i\"" "\"$FS\"" "\"$j\"" "\"$NW\"" "\"$K\"" "\"$PVAL\"" "\"$3\"" "\"$4\""
  done
done
