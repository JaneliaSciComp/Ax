using MAT
#using Debug
#using Profile
using Distributions
using NumericExtensions

#@iprofile begin
function ftest(data::Array{Float32,2}, tapers::Array{Float32,2}, p::Float64,
    fft_in::Array{Float32,3}, fft_out::Array{Complex{Float32},3}, fft_plan::Any)

  (NC,C)=size(data)
  (NK,K)=size(tapers);
  N=NC

  #f = zeros(N)  # not sure above is correct

  Kodd = [1:2:K]
  Keven = [2:2:K]

  tapers2 = broadcast(*, tapers, ones(Float32,1,1,1))  # f x K x 1
  data = broadcast(*, data, ones(Float32,1,1,1))       # f x C x 1
  data = permutedims(data,(1,3,2))                     # f x 1 x C
  broadcast!(*, fft_in, data, tapers2)                 # f x K x C
  FFTW.execute(fft_plan.plan, fft_in, fft_out)
  J = fft_out                                    # f x K x C

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

  #sig=invlogcdf(FDist(2,2*K-2),log(1-p/N))
  #var=den./(K*H0sq)
  #sd=sqrt(var)

  #A=A*Fs

  #return Fval, A, f, sig, sd
  return Fval

end
#end #profile

#@iprofile begin
function brown_puckette(x,k,fs, fft_in2, fft_out2, fft_plan2)
  nfft = length(x)
  #X=np.fft.fft(x)
  fft_in2 = x;
  FFTW.execute(fft_plan2.plan, fft_in2, fft_out2)
  X = fft_out2
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

function do_it(params)
  DIROUT = params[1]
  FILEINs = params[2]
  t_now = params[3]
  NW = params[4]
  K = params[5]
  PVAL = params[6]
  FS = params[7]
  NFFT = params[8]
  CHUNK = params[9]
  tapers = params[10]
  offset = params[11]
  offset2 = params[12]

  NCHANNELS = length(FILEINs)
  sig=invlogcdf(FDist(2,2*K-2),log(1-PVAL/NFFT))

  fft_in = zeros(Float32, NFFT, K, NCHANNELS)
  fft_out = zeros(Complex64, NFFT>>1+1, K, NCHANNELS)
  fft_plan = FFTW.Plan(fft_in, fft_out, 1, FFTW.PATIENT, FFTW.NO_TIMELIMIT)
  fft_in2 = zeros(Float32, NFFT)
  fft_out2 = zeros(Complex64, NFFT>>1+1)
  fft_plan2 = FFTW.Plan(fft_in2, fft_out2, 1, FFTW.PATIENT, FFTW.NO_TIMELIMIT)

  #dd[1:(NFFT>>1),:]=dd[(end-NFFT>>1+1):end,:]
  NSAMPLES = (NFFT>>1)*(CHUNK+1);
  dd = Array(Float32, NSAMPLES, NCHANNELS)
  for i = 1:NCHANNELS
    fid = open(string(DIROUT,"/",FILEINs[i]),"r")
    seek(fid, ((t_now+offset*CHUNK)*(NFFT>>1)+offset2)*4)
    tmp = read(fid, Float32, NSAMPLES)
    if (length(tmp) < NSAMPLES)
      tmp = [tmp, zeros(Float32, NSAMPLES-length(tmp))]
    end
    dd[:,i]=tmp
    #dd[(NFFT>>1+1):end,i] = tmp
  end

  idx={}
  for j = 1:CHUNK
    ddd=dd[[1:NFFT]+(j-1)*(NFFT>>1),:]
    #(F, A, f, sig, sd) = ftest(ddd, tapers, FS, PVAL, fft_in, fft_out, fft_plan)
    F = ftest(ddd, tapers, PVAL, fft_in, fft_out, fft_plan)
    for l = 1:NCHANNELS
      tmp = 1 + find(F[2:end-1,l] .> sig)
      for m = 1:length(tmp)
        freq, amp = brown_puckette(ddd[:,l],tmp[m],FS, fft_in2, fft_out2, fft_plan2)
        push!(idx, [float64(j)+offset*CHUNK, freq, amp, float64(l)])
      end
    end
  end
  return idx
end
