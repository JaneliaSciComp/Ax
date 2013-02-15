#!/bin/bash

# cd to folder containing cluster.sh and then execute:
#   ./cluster.sh  full_path_to_parameters  full_path_to_data [start(s) stop(s)]
# data can either be a folder of sessions or a single session's base filename
# for example
#   ./cluster.sh  ./ultrasonic_parameters.m  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/
#   ./cluster.sh  ./ultrasonic_parameters.m  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/Test_B_1

# FS, NW, K, PVAL, and NFFT specify the parameters to each call of ax().
# if any are a scalar the same value is used for each call.  those which are
# arrays must have the same length.

if [ $# -ne 2 ] && [ $# -ne 4 ]
then
  echo invalid arguments
  exit
fi

eval $( sed "s/\[/\\(/g" $1 | sed "s/\]/\\)/g" )

# get maximum length of params
tmp=(${#FS[@]} ${#NW[@]} ${#K[@]} ${#PVAL[@]} ${#NFFT[@]})
tmp=$(echo ${tmp[@]} | awk -v RS=" " '1' | sort -nr | head -1)

# expand scalars to vector of above max len
if [ ${#FS[@]} -lt $tmp ] ; then
  for i in $(seq 1 $(($tmp - 1)))
  do FS[i]=${FS[0]}; done
fi

if [ ${#NW[@]} -lt $tmp ] ; then
  for i in $(seq 1 $(($tmp - 1)))
  do NW[i]=${NW[0]}; done
fi

if [ ${#K[@]} -lt $tmp ] ; then
  for i in $(seq 1 $(($tmp - 1)))
  do K[i]=${K[0]}; done
fi

if [ ${#PVAL[@]} -lt $tmp ] ; then
  for i in $(seq 1 $(($tmp - 1)))
  do PVAL[i]=${PVAL[0]}; done
fi

if [ ${#NFFT[@]} -lt $tmp ] ; then
  for i in $(seq 1 $(($tmp - 1)))
  do NFFT[i]=${NFFT[0]}; done
fi

#echo ${FS[@]}
#echo ${NW[@]}
#echo ${K[@]}
#echo ${PVAL[@]}
#echo ${NFFT[@]}

# launch one instance of ax() per set of params
if [ -d $2 ] ; then
  ii=$(ls -1 $2/*.ch* | sed s/.ch[0-9]// | uniq)
  dir_name=$2
else
  ii=$2
  dir_name=$(dirname $2)
fi
for i in $ii
do
#  echo $i
  job_name=$(basename $i)
  for j in $(seq 0 $((${#NFFT[@]} - 1)))
  do
#    echo $j
    qsub -N "$job_name-$j" -pe batch 8 -b y -j y -cwd -o "$dir_name/$job_name-$j.log" -V ./cluster2.sh "\"${FS[j]}\"" "\"${NFFT[j]}\"" "\"${NW[j]}\"" "\"${K[j]}\"" "\"${PVAL[j]}\"" "\"$i\"" "\"$j\"" "\"$3\"" "\"$4\""
    sleep 1
  done
done
