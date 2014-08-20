#using MAT
using WAV
#using Debug
#using Profile
using NumericExtensions
using HDF5

require("lock.jl")

#@iprofile begin
function ftest(data::Array{Float32,2}, tapers::Array{Float32,2}, p::Float64)

  (NC,C)=size(data)
  (NK,K)=size(tapers);
  N=NC

  Kodd = [1:2:K]
  Keven = [2:2:K]

  tapers2 = broadcast(*, tapers, ones(Float32,1,1,1))  # f x K x 1
  data = broadcast(*, data, ones(Float32,1,1,1))       # f x C x 1
  data = permutedims(data,(1,3,2))                     # f x 1 x C
  data = broadcast(*, data, tapers2)                   # f x K x C

  # use slow rfft for now b/c fast plan_rfft is buggy and in flux
  J = rfft(data,1)                                     # f x K x C

  Jp = J[:,Kodd,:]                               # f x K x C
  H0 = sum(tapers2[:,Kodd],1)                    # 1 x K
  H0sq = sumsq(H0)                               # 1
  JpH0 = squeeze(sum(broadcast(*, Jp, H0),2),2)  # f x C
  A = JpH0./H0sq                                 # f x C
  Kp = size(Jp,2)
  Ap = broadcast(*, A, ones(Float32,1,1,1))      # f x C x 1
  Ap = permutedims(Ap,(1,3,2))                   # f x 1 x C
  Jhat = broadcast(*, Ap, H0)                    # f x K x C
  num=(K-1)*abs2(A).*H0sq
  den=squeeze(sum(sqrdiff(Jp,Jhat),2),2)+squeeze(sum(abs2(J[:,Keven,:]),2),2)
  Fval=num./den

  return Fval

end
#end #profile

#@iprofile begin
function brown_puckette(x,k,fs)
  nfft = length(x)
  X=rfft(x,1)
  Xh0 = 0.5*(X[k]-0.5*X[k+1]-0.5*X[k-1])
  Xh1 = 0.5*exp(im*2*pi*(k-1)/nfft) *
      (X[k] - 0.5*exp(im*2*pi/nfft)*X[k+1] - 0.5*exp(-1im*2*pi/nfft)*X[k-1])
  phi0 = atan2(imag(Xh0),real(Xh0))
  phi1 = atan2(imag(Xh1),real(Xh1))
  if (phi1-phi0)<0
    phi1 = phi1+2*pi
  end
  freq = (phi1-phi0)*fs/(2*pi)

  period = fs/freq
  last = int(floor(period * floor(length(x)/period)))
  if last>0
    real_part = mean(x[1:last] .* cos([1:last].*(2*pi/period)))
    imag_part = mean(x[1:last] .* sin([1:last].*(2*pi/period)))
    amp = 2*abs(real_part + im*imag_part)
  else
    amp = NaN
  end
  return freq,amp
end
#end #profile

#@debug function do_it(params)
function do_it(params)
  FILEPATH = params[1]
  FILENAMES = params[2]
  FILENAME = params[3]
  FILETYPE = params[4]
  FILELEN_TIC = params[5]
  NCHANNELS = params[6]
  t = params[7]
  NWINDOWS_PER_WORKER = params[8]
  NFFT = params[9]
  PVAL = params[10]
  FS = params[11]
  tapers = params[12]
  sig = params[13]

  NSAMPLES = (NFFT>>1)*(NWINDOWS_PER_WORKER+1);
  if FILETYPE=="ch"
    for i = 1:NCHANNELS
      filei=string(FILEPATH,"/",FILENAMES[i])
      fid = open(filei,"r")
      seek(fid, t*4)
      tmp2 = read(fid, Float32, min(NSAMPLES, FILELEN_TIC-t))
      if(i==1)  dd=similar(tmp2, (size(tmp2,1), NCHANNELS))  end
      dd[:,i]=tmp2;
      close(fid)
    end
  else #if FILETYPE=="wav"
    dd, fs, nbits, extra =
        wavread("$FILENAME.wav"; subrange=(t+1):min(t+NSAMPLES, FILELEN_TIC))
    dd=float32(dd)
  end
  if (size(dd,1) < NSAMPLES)
    dd = [dd, zeros(Float32, (NSAMPLES-size(dd,1),NCHANNELS))]
  end

  idx={}
  for j = 1:NWINDOWS_PER_WORKER
    ddd=dd[[1:NFFT].+(j-1)*(NFFT>>1),:]
    F = ftest(ddd, tapers, PVAL)
#@bp
    for l = 1:NCHANNELS
      tmp = 1 .+ find(F[2:end-1,l] .> sig)
      for m = 1:length(tmp)
        freq, amp = brown_puckette(ddd[:,l],tmp[m],FS)
        push!(idx, [t/(NFFT>>1)+float64(j), freq, amp, float64(l)])
      end
    end
  end
  @printf("done with %.1f sec\n",t/FS);
  return [x[y] for x in idx, y=1:4]
end
