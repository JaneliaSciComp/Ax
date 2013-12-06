#!/home/arthurb/bin/anaconda/bin/python

# python ax1.py params_file FILEIN FILEOUT
# python ax1.py params_file FILEIN FILEOUT START STOP
# python ax1.py FS NFFT NW K PVAL FILEIN FILEOUT
# python ax1.py FS NFFT NW K PVAL FILEIN FILEOUT START STOP
#
# analyze a set of time series with multi-taper spectral analysis and
# create a sparse matrix of just the time-frequency pixels whose F-test
# passes PVAL.
#
# typical usage consists of one or more input files being analyzed by one
# or more parameter sets.  for example, four microphone recordings of the
# same vocalizing mouse analyzed with three different NFFTs and the same
# NW, K, and PVAL.  <filename>.ch[1-4] yield <filename>-[1-3].ax
#
# FS: sampling rate in Hertz
# NFFT: FFT window size in seconds, rounds up to the next power of 2 tics
# NW: multi-taper time-bandwidth product
# K: number of tapers
# PVAL: F-test p-val threshold
# FILEIN: the base filename and path of [0-9].wav files with a single channel each,
#     or .ch[0-9] files containing float32s
# FILEOUT: an integer to append to FILEIN to differentiate parameter sets used
# START,STOP: optional time range, in seconds
#
# output is a binary file with a time x frequency x amplitude x channel
#     array of hot pixels
#
# python ax1.py 'ultrasonic_params.txt' 'urine' '1'
# python ax1.py 200e3 0.001 15 29 0.01 'urine' '1'
# python ax1.py 450450 0.001 15 29 0.01 0 30 'groundtruth' '1'

# /home/arthurb/bin/anaconda/bin/kernprof.py -l -v ax1.py 450450 0.00025 22 43 0.01 /groups/egnor/egnorlab/ben/Test_D_1 7 0 4

from ax1b import do_it, nextpow2
import struct 
import time
import numpy as np
import glob
import sys
import os
from multiprocessing import Pool, cpu_count
#import pdb
import math
from scipy import stats
import pyfftw
from dpss import dpss
import wave

if __name__ == "__main__":

  if (len(sys.argv)!=4) and (len(sys.argv)!=6) and (len(sys.argv)!=8) and (len(sys.argv)!=10):
    print('invalid args')
    exit


  tstart=time.time()

  if (len(sys.argv)<8):
    execfile(sys.argv[1])
    FILEIN=sys.argv[2]
    FILEOUT=sys.argv[3]
  else:
    FS=sys.argv[1]
    NFFT=sys.argv[2]
    NW=sys.argv[3]
    K=sys.argv[4]
    PVAL=sys.argv[5]
    FILEIN=sys.argv[6]
    FILEOUT=sys.argv[7]
  if ((len(sys.argv)==6) or (len(sys.argv)==10)):
    START=sys.argv[-2]
    STOP=sys.argv[-1]

  if (isinstance(FS,str)):
    FS = int(FS)
  if (isinstance(NFFT,str)):
    NFFT = float(NFFT)
  if (isinstance(NW,str)):
    NW = int(NW)
  if (isinstance(K,str)):
    K = int(K)
  if (isinstance(PVAL,str)):
    PVAL = float(PVAL)
  if ((len(sys.argv)==6) or (len(sys.argv)==10)):
    if (isinstance(START,str)):
      START = float(START)
    if (isinstance(STOP,str)):
      STOP = float(STOP)

  VERSION=1

  SUBSAMPLE=1
  NWORKERS=cpu_count()

  FS=int(FS/SUBSAMPLE);

  NFFT=int(nextpow2(NFFT*FS))  # convert to ticks

  NWINDOWS_PER_WORKER=int(12*256*1000/NFFT)  # NFFT/2 ticks

  FIRST_MT=float('nan')
  LAST_MT=float('nan')
  FRACTION_MT=float('nan')

  tapers,eig = dpss(NFFT, NW, K)
  tapers = np.array(tapers, dtype=np.float32)
  #tapers = tapers * np.sqrt(FS)

  f=np.array(range(0,NFFT//2+1))*FS/NFFT
  df=f[1]-f[0];

  DIROUT=os.path.dirname(FILEIN);
  FILEINs=sorted(glob.glob(FILEIN+'.ch*'));
  FILE_TYPE=1
  if (len(FILEINs)==0):
    FILEINs=sorted(glob.glob(FILEIN+'*.wav'));
    FILE_TYPE=2
    if (len(FILEINs)==0):
      print(["can't find any .wav or .ch files with basename '"+FILEIN]);
      exit
  NCHANNELS=len(FILEINs);

  REMAP=list();
  for i in range(0,NCHANNELS):
    filei=os.path.join(DIROUT,FILEINs[i])
    if FILE_TYPE==1:
      try:
        fid=open(filei,'rb')
      except:
        print(["can't open file '"+filei+"'"])
        exit
      fid.seek(0,2);
      FILE_LEN=fid.tell()/4/FS;
      fid.close()
      REMAP.append(FILEINs[i][-1]);
    if FILE_TYPE==2:
      try:
        fid=wave.open(filei,'rb')
      except:
        print(["can't open file '"+filei+"'"])
        exit
      FILE_LEN=fid.getnframes()/FS
      fid.close();
      REMAP.append(FILEINs[i][-5]);
    if 'START' not in locals():
      tmp=FILE_LEN*FS/(NFFT//2)-1
      print('Processing {:.3g} min = {:.3g} windows = {:3g} chunks of data in {:s}'.format(FILE_LEN/60, tmp, tmp/NWINDOWS_PER_WORKER, FILEINs[i]));
      t_offset_tic=0;
      t_now_sec=0;
    else:
      tmp=(STOP-START)*FS/(NFFT//2)-1
      print('Processing {:.3g} min = {:.3g} windows = {:3g} chunks of data in {:s}'.format((STOP-START)/60, tmp, tmp/NWINDOWS_PER_WORKER, FILEINs[i]));
      t_offset_tic=round(START*FS);
      t_now_sec=START;

  fid_out=open(FILEIN+'-'+FILEOUT+'.ax','wb')
  # L=8 bytes on 64-bit systems
  fid_out.write(struct.pack('B',VERSION))
  fid_out.write(struct.pack('B',SUBSAMPLE))
  fid_out.write(struct.pack('B',0))
  fid_out.write(struct.pack('I',FS))
  fid_out.write(struct.pack('I',NFFT))
  fid_out.write(struct.pack('H',NW))
  fid_out.write(struct.pack('H',K))
  fid_out.write(struct.pack('d',PVAL))
  fid_out.write(struct.pack('d',df))

  t_now=0
  tloop=time.time()
  pool=Pool()

  while ((t_now_sec<FILE_LEN) and (('STOP' not in locals()) or (t_now_sec<STOP))):
    if ((time.time()-tloop)>10):
      tmp=t_now_sec
      tmp2=0
      if 'START' in locals():
        tmp=tmp-START
        tmp2=START
      if 'STOP' in locals():
        tmp=tmp/(STOP-tmp2)
      else:
        tmp=tmp/(FILE_LEN-tmp2)
      print('{:d} sec processed;  {:d}% done'.format(int(round(t_now_sec-tmp2)),int(round(100*tmp))))
      tloop=time.time()

    #idx=map(do_it, \
    idx=pool.map(do_it, \
       [(DIROUT, FILEINs, t_now, NW,K,PVAL,FS,NFFT, NWINDOWS_PER_WORKER, tapers, x, t_offset_tic, FILE_TYPE, round(FILE_LEN*FS)) for x in range(0,NWORKERS)])
    for i in idx:
      for j in i:
        fid_out.write(struct.pack('dddd', \
            float(t_now)+j[0], j[1], j[2], float(REMAP[j[3]])))

    t_now_sec = t_now_sec+float(NFFT//2)/FS*NWORKERS*NWINDOWS_PER_WORKER
    t_now = t_now+NWORKERS*NWINDOWS_PER_WORKER
 
  fid_out.write('Z'.encode('ascii'))
  fid_out.close()
 
  tstop = time.time() - tstart
  print('Run time was {:.3g} minutes.'.format(tstop/60))

  pool.close()
