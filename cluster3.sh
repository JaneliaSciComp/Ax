#!/bin/bash

date

hostname

export MCR_CACHE_VERBOSE=1

if [ -d /scratch/$USER ] ; then
  export MCR_CACHE_ROOT=/scratch/$USER/mcr_cache_root.$JOB_ID
else
  export MCR_CACHE_ROOT=~/mcr_cache_root.$JOB_ID
fi

if [ -d MCR_CACHE_ROOT ] ; then
  echo Deleting pre-existing MCR_CACHE_ROOT
  rm -rf $MCR_CACHE_ROOT
fi

mkdir $MCR_CACHE_ROOT

run_ax2.sh /usr/local/matlab-2014a $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}

rm -rf $MCR_CACHE_ROOT

date
