#!/bin/bash

# cd to folder containing cluster.sh and then execute:
#   ./cluster.sh  <full_path_to_parameters>  <full_path_to_data> <multitaper_flag> [<start(s)> <stop(s)>]
#
# data can either be a folder of sessions or a single session's base filename
# use 1 for multitaper_flag to run the fourier analysis in ax1 and the heuristics in ax2
# use 0 for multitaper_flat to run just ax2 on pre-exiting .ax files
#
# for example
#   ./cluster.sh  ./ultrasonic_parameters.m  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/ 0
#   ./cluster.sh  ./ultrasonic_parameters.m  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/Test_B_1 1

# FS, NW, K, PVAL, and NFFT specify the parameters to each call of ax1().
# if any are a scalar the same value is used for each call.  those which are
# arrays must have the same length.

if [ $# -ne 3 ] && [ $# -ne 5 ]
then
  echo invalid arguments
  exit
fi

if [ $3 -eq 1 ]
then
  echo starting ax1

  eval $( sed "s/\[/\\(/g" $1 | sed "s/\]/\\)/g" | dos2unix )

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

  if [ -d $2 ] ; then
    ii=$(ls -1 $2/*.ch* | sed s/.ch[0-9]// | uniq)
    dir_name=$2
  else
    ii=$2
    dir_name=$(dirname $2)
  fi
  k=0
  for i in $ii
  do
  #  echo $i
    job_name=$(basename $i)
    for j in $(seq 0 $((${#NFFT[@]} - 1)))
    do
      k=$((k+1))
  #    echo $j
      qsub -N "$job_name-$j" \
          -pe batch 8 \
          -b y -j y -cwd -o "$dir_name/$job_name-$j.log" \
          -V ./cluster2.sh "\"${FS[j]}\"" "\"${NFFT[j]}\"" "\"${NW[j]}\"" "\"${K[j]}\"" "\"${PVAL[j]}\"" "\"$i\"" "\"$j\"" "\"$4\"" "\"$5\""
#          -pe batch 12 \
#          -l r620=true \
    done
  done

  echo waiting for ax1 to finish
  sleep 1m
  while [ $(grep "Run time was " $2*.log | wc -l) -ne $k ]
  do
    sleep 1m
  done

fi

echo starting ax2
if [ -d $2 ] ; then
  ii=$(ls -1 $2/*.ch* | sed s/.ch[0-9]// | uniq)
  dir_name=$2
else
  ii=$2
  dir_name=$(dirname $2)
fi
for i in $ii
do
  job_name=$(basename $i)
  qsub -N "$job_name" \
      -pe batch 8 \
      -b y -j y -cwd -o "$dir_name/$job_name.log" \
      -V ./cluster3.sh "\"$1\"" "\"$i\""
#      -pe batch 12 \
#      -l r620=true \
done
echo check the .log file to see when ax2 finishes
echo goodbye
