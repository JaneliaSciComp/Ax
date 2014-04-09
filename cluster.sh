#!/bin/bash

# cd to folder containing cluster.sh and then execute:
#   ./cluster.sh  <full_path_to_parameters>  <full_path_to_data> <multitaper_flag> [<start(s)> <stop(s)>]
#
# data can either be a folder of sessions or a single session's base filename
# for multitaper_flag,
#   use 1 to run just the fourier analysis in ax1
#   use 2 to run just ax2 on pre-exiting .ax files
#   use 3 to run both the fourier analysis in ax1 and the heuristics in ax2
#
# for example
#   ./cluster.sh  ./ultrasonic_parameters.m  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/ 0
#   ./cluster.sh  ./ultrasonic_parameters.m  /groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/Test_B_1 1

# FS, NW, K, PVAL, and NFFT specify the parameters to each call of ax1().
# if any are a scalar the same value is used for each call.  those which are
# arrays must have the same length.

if [ $# -ne 3 ] && [ $# -ne 5 ] ; then
  echo invalid arguments
  exit
fi

if [ -d $2 ] ; then
  ii=$(ls -1 ${2%/}/*.ch* | sed s/\.ch[0-9]*// | uniq)
  dir_name=$2
else
  ii=$2
  dir_name=$(dirname $2)
fi
#echo $ii

if [ $(($3 & 1)) -gt 0 ] ; then
  echo starting ax1

  eval $( sed "s/\[/\\(/g" $1 | sed "s/\]/\\)/g" | dos2unix )

  # get maximum length of params
  tmp=(${#FS[@]} ${#NW[@]} ${#K[@]} ${#PVAL[@]} ${#NFFT[@]})
  tmp=$(echo ${tmp[@]} | awk -v RS=" " '1' | sort -nr | head -1)

  # expand scalars to vector of above max len
  if [ ${#FS[@]} -lt $tmp ] ; then
    for i in $(seq 1 $(($tmp - 1))) ; do
      FS[i]=${FS[0]}
    done
  fi

  if [ ${#NW[@]} -lt $tmp ] ; then
    for i in $(seq 1 $(($tmp - 1))) ; do
      NW[i]=${NW[0]}
    done
  fi

  if [ ${#K[@]} -lt $tmp ] ; then
    for i in $(seq 1 $(($tmp - 1))) ; do
      K[i]=${K[0]}
    done
  fi

  if [ ${#PVAL[@]} -lt $tmp ] ; then
    for i in $(seq 1 $(($tmp - 1))) ; do
      PVAL[i]=${PVAL[0]}
    done
  fi

  if [ ${#NFFT[@]} -lt $tmp ] ; then
    for i in $(seq 1 $(($tmp - 1))) ; do
      NFFT[i]=${NFFT[0]}
    done
  fi

  #echo ${FS[@]}
  #echo ${NW[@]}
  #echo ${K[@]}
  #echo ${PVAL[@]}
  #echo ${NFFT[@]}

  k=0
  for i in $ii ; do
    #echo $i
    job_name=$(basename $i)
    for j in $(seq 0 $((${#NFFT[@]} - 1))) ; do
      k=$((k+1))
      #echo $j
      qsub -N "$job_name-$j" \
          -pe batch 16 \
          -b y -j y -cwd -o "$dir_name/$job_name-$j.log" \
          -V cluster2.sh "\"${FS[j]}\"" "\"${NFFT[j]}\"" "\"${NW[j]}\"" "\"${K[j]}\"" "\"${PVAL[j]}\"" "\"$i\"" "\"$j\"" "\"$4\"" "\"$5\""
    done
  done
fi

if [ $3 -eq 3 ] ; then
  echo waiting for ax1 to finish
  # this won't work if .log files already exist from a previous run
  iii=$(sed s/$/*.log/ <<< "$ii")
  #echo $iii
  sleep 1m
  while [ $(grep "Run time was " $iii | wc -l) -ne $k ] ; do
    sleep 1m
  done
fi

if [ $(($3 & 2)) -gt 0 ] ; then
  echo starting ax2
  for i in $ii ; do
    job_name=$(basename $i)
    qsub -N "$job_name" \
        -b y -j y -cwd -o "$dir_name/$job_name.log" \
        -V cluster3.sh "\"$1\"" "\"$i\""
  done
  echo check the .log file to see when ax2 finishes
  echo goodbye
fi
