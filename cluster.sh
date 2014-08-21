#!/bin/bash

# cd to folder containing cluster.sh and then execute:
#   ./cluster.sh <full_path_to_parameters> <full_path_to_data> <multitaper_flag> [<start(s)> <stop(s)>]
#
# <full_path_to_data> can either be a folder of .wav or .ch files, or a single .wav file,
#   or a group of .ch files w/o the suffix
#
# for multitaper_flag,
#   use 1 to run just the fourier analysis in ax1
#   use 2 to run just the heuristics in ax2 on pre-exiting ax1 output
#   use 3 to run both the fourier analysis in ax1 and the heuristics in ax2
#
# for example
#   ./cluster.sh ./parameters.txt /lab/postdoc2/expA/data 3
#   ./cluster.sh ./parameters.txt /lab/gradstudent3/expC/trialdata/foo.wav 1 10 30
#
# FS, NW, K, PVAL, and NFFT specify the parameters to each call of ax1.
# if any are a scalar the same value is used for each call.  those which are
# arrays must have the same length.

if [ $# -ne 3 ] && [ $# -ne 5 ] ; then
  echo invalid arguments
  exit
fi

# parse parameters file
eval $( sed "s/\[/\\(/g" $1 | sed "s/\]/\\)/g" | dos2unix )

# get list of data files
if [ -d $2 ] ; then
  ii=$(ls -1 ${2%/}/*.wav ${2%/}/*.WAV 2> /dev/null)
  ii="${ii} $(ls -1 ${2%/}/*.ch* 2> /dev/null | sed s/\.ch[0-9]$// | uniq)"
  dir_name=$2
else
  ii=$2
  dir_name=$(dirname $2)
fi
#echo $ii

# run ax1
if [ $(($3 & 1)) -gt 0 ] ; then
  echo starting ax1

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
  logfiles=""
  for i in $ii ; do
    #echo $i
    i_wavless=${i%.wav}
    i_wavless=${i_wavless%.WAV}
    job_name=$(basename $i_wavless)
    for j in $(seq 0 $((${#NFFT[@]} - 1))) ; do
      #echo $j
      k=$((k+1))
      logfiles="${logfiles} ${i_wavless}-${j}.log"
      qsub -N "ax$job_name-$j" \
          -pe batch 16 \
          -b y -j y -cwd -o "$dir_name/$job_name-$j.log" \
          -V cluster1.sh "\"${FS[j]}\"" "\"${NFFT[j]}\"" "\"${NW[j]}\"" "\"${K[j]}\"" "\"${PVAL[j]}\"" "\"$i\"" "\"$j\"" "\"$4\"" "\"$5\""
    done
  done
fi

if [ $3 -eq 3 ] ; then
  echo waiting for queue to flush
  while ! (ls ${logfiles} > /dev/null 2>&1) ; do
    sleep 1m
  done
  echo waiting for ax1 to finish
  # this won't work if .log files already exist from a previous run
  while [ $(grep "Run time was " ${logfiles} | wc -l) -ne $k ] ; do
    sleep 1m
  done
fi

# run ax2
if [ $(($3 & 2)) -gt 0 ] ; then
  echo starting ax2
  for i in $ii ; do
    i_wavless=${i%.wav}
    i_wavless=${i_wavless%.WAV}
    job_name=$(basename $i_wavless)
    axfiles=""
    for j in $(seq 0 $((${#NFFT[@]} - 1))) ; do
      axfiles="${axfiles} ${i_wavless}-${j}.ax"
    done
    qsub -N "ax$job_name" \
        -b y -j y -cwd -o "$dir_name/$job_name.log" \
        -V cluster2.sh "\"${frequency_low}\"" "\"${frequency_high}\"" "\"${convolution_size[*]}\"" \
        "\"${minimum_object_area}\"" "\"${merge_harmonics}\"" "\"${merge_harmonics_overlap}\"" \
        "\"${merge_harmonics_ratio}\"" "\"${merge_harmonics_fraction}\"" \
        "\"${minimum_vocalization_length}\"" "\"${channels[*]}\"" "\"${axfiles}\"" "\"${i_wavless}\""
  done
  echo check the .log file to see when ax2 finishes
  echo goodbye
fi
