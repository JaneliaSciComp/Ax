#from numba import autojit
import pyfftw
import struct 
import numpy as np
import glob
import sys
import os
#import pdb
import math
from scipy import stats
from dpss import dpss

def nextpow2(i):
    n = 2
    while n < i: n = n * 2
    return n

#@autojit
#@profile
def ftest(data, tapers, p, fft_in, fft_out, fft_plan):
  tmp = data.shape
  N = tmp[0]
  C = tmp[1]
  tmp = tapers.shape
  K = tmp[1]

  #df=float(Fs)/N
  #f=np.arange(0,Fs,df)
  #f=f[0:(N//2+1)]

  Kodd = np.arange(0,K,2)
  Keven = np.arange(1,K,2)

  tapers2 = tapers[:,:,np.newaxis]
  data = data[:,np.newaxis,:]
  fft_in[:] = data * tapers2
  fft_plan()
  J=fft_out

  Jp = J[0:(N//2+1),Kodd,:]         # f x K x C
  H0 = np.sum(tapers[:,Kodd],0)     # K
  H0sq = np.empty((N//2+1, C))      # f x C
  H0sq.fill(sum(H0*H0))             # much slower to not fill it, strangely
  H0 = H0[np.newaxis,:,np.newaxis]  # 1 x K x 1
  JpH0 = np.sum(Jp*H0 ,1)           # f x C
  A = JpH0/H0sq                     # f x C
  tmp = Jp.shape
  Kp = tmp[1]
  Ap = A[:,:,np.newaxis]            # f x C x 1
  Ap = np.transpose(Ap,(0,2,1))     # f x 1 x C
  Jhat = Ap * H0;                   # f x K x C

  num=(K-1)*(A.real*A.real+A.imag*A.imag)*H0sq
  den1 = Jp-Jhat
  den1 = den1.real*den1.real + den1.imag*den1.imag
  den1 = np.sum(den1,1)
  den2 = J[0:(N//2+1),Keven,:]
  den2 = den2.real*den2.real + den2.imag*den2.imag
  den2 = np.sum(den2,1)
  den = den1 + den2
  Fval=num/den

  #sig=stats.f.ppf((1-p/N),2,2*K-2)
  #var=den/(K*H0sq)
  #sd=np.sqrt(var)

  #A=A*Fs

  #return Fval, A, f, sig, sd
  return Fval


#% from charpentier (1986) and brown and puckette (1993; JASA)
def brown_puckette(x,k,fs, fft_in2, fft_out2, fft_plan2):
  nfft=len(x)
  fft_in2[:]=x
  fft_plan2()
  X=fft_out2
  Xh0=0.5*(X[k]-0.5*X[k+1]-0.5*X[k-1])
  Xh1=0.5*np.exp(1j*2*np.pi*(k-0)/nfft)* \
      (X[k] - 0.5*np.exp(1j*2*np.pi/nfft)*X[k+1] \
            - 0.5*np.exp(-1j*2*np.pi/nfft)*X[k-1])
  phi0=math.atan2(np.imag(Xh0),np.real(Xh0))
  phi1=math.atan2(np.imag(Xh1),np.real(Xh1))
  if((phi1-phi0)<0):
    phi1=phi1+2*np.pi
  freq=(phi1-phi0)*fs/(2*np.pi)

  amp = np.nan
  if freq>0:
    period = fs/freq
    last = int(np.floor(period * np.floor(len(x)/period)))
    if last>0:
      real_part = np.mean(np.multiply(x[0:last], \
          np.cos(np.array(range(0,last), dtype=np.float32)*(2*np.pi/period))))
      imag_part = np.mean(np.multiply(x[0:last], \
          np.sin(np.array(range(0,last), dtype=np.float32)*(2*np.pi/period))))
      amp = 2*np.abs(complex(real_part, imag_part))

  return freq,amp

def do_it(params):
  DIROUT = params[0]
  FILEINs = params[1]
  t_now = params[2]
  NW = params[3]
  K = params[4]
  PVAL = params[5]
  FS = params[6]
  NFFT = params[7]
  CHUNK = params[8]
  tapers = params[9]
  offset = params[10]
  offset2 = params[11]

  NCHANNELS = len(FILEINs)
  sig=stats.f.ppf((1-PVAL/NFFT),2,2*K-2)

  fft_in = pyfftw.n_byte_align_empty((NFFT,K,NCHANNELS), 16, 'float32')
  fft_out = pyfftw.n_byte_align_empty((NFFT//2+1,K,NCHANNELS), 16, 'complex64')
  fft_plan = pyfftw.FFTW(fft_in, fft_out, axes=(0,), flags=('FFTW_PATIENT',))
  fft_in2 = pyfftw.n_byte_align_empty(NFFT, 16, 'float32')
  fft_out2 = pyfftw.n_byte_align_empty(NFFT//2+1, 16, 'complex64')
  fft_plan2 = pyfftw.FFTW(fft_in2, fft_out2, flags=('FFTW_PATIENT',))

  NSAMPLES = NFFT//2*(CHUNK+1)
  dd = np.empty([NSAMPLES, NCHANNELS], dtype=np.float32)
  for i in range(0,NCHANNELS):
    fid=open(os.path.join(DIROUT,FILEINs[i]),'rb')
    fid.seek(((t_now+offset*CHUNK)*(NFFT//2)+offset2)*4,0);
    tmp = fid.read(NSAMPLES*4)
    tmp = np.array(struct.unpack(str(int(len(tmp)/4))+'f', tmp), dtype=np.float32)
    if (len(tmp) < NSAMPLES):
      tmp = np.concatenate((tmp, np.zeros(NSAMPLES-len(tmp), dtype=np.float32)))
    dd[:,i] = tmp
    fid.close()

  idx=list()
  for j in range(0,CHUNK):
    ddd=dd[np.array(range(0,NFFT))+j*NFFT//2,:]
    #F, A, f, sig, sd = ftest(ddd, tapers, FS, PVAL, fft_in, fft_out, fft_plan)
    F = ftest(ddd, tapers, PVAL, fft_in, fft_out, fft_plan)
    for l in range(0,NCHANNELS):
      tmp=[i+1 for (i,v) in enumerate(F[1:-1,l]) if v>sig]
      for m in range(0,len(tmp)):
        freq,amp = brown_puckette(ddd[:,l],tmp[m],FS, fft_in2, fft_out2, fft_plan2)
        idx.append((j+offset*CHUNK, freq, amp, l))
  return idx
