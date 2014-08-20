#!/home/arthurb/src/julia/julia

# julia -p NWORKERS ax1.jl params_file FILENAME SUFFIX
# julia -p NWORKERS ax1.jl params_file FILENAME SUFFIX START STOP
# julia -p NWORKERS ax1.jl FS NFFT NW K PVAL FILENAME SUFFIX
# julia -p NWORKERS ax1.jl FS NFFT NW K PVAL FILENAME SUFFIX START STOP
#
# analyze a set of time series with multi-taper spectral analysis and
# create a sparse matrix of just the time-frequency pixels whose F-test
# passes PVAL.
#
# typical usage consists of one or more input files being analyzed by one
# or more parameter sets.  for example, four microphone recordings of the
# same vocalizing mouse analyzed with three different NFFTs and the same
# NW, K, and PVAL.  <filename>.wav yield <filename>-[1-3].ax
#
# FS: sampling rate in Hertz
# NFFT: FFT window size in seconds, rounds up to the next power of 2 tics
# NW: multi-taper time-bandwidth product
# K: number of tapers
# PVAL: F-test p-val threshold
# FILENAME: the full path to a single .wav file containing all channels,
#     or of .ch[0-9] files each with a single channel of float32s
# SUFFIX: a string to append to FILENAME to differentiate parameter sets used
# START,STOP: optional time range, in seconds
#
# output is a binary file with a time x frequency x amplitude x channel
#     array of hot pixels
#
# julia -p 12 ax1.jl 'ultrasonic_params' 'urine' '1'
# julia -p 12 ax1.jl 200e3 32 15 29 0.01 'urine' '1'
# julia -p 12 ax1.jl 450450 32 15 29 0.01 0 30 'groundtruth' '1'

#using MAT
using DSP
using WAV
using HDF5
#using Debug
#using Profile
using Distributions

#using ClusterManagers
#addprocs(10, cman=SGEManager())
#require("/home/arthurb/projects/egnor/ax/ax1b.jl")  # full path on cmd line too

require("ax1b.jl")

#@debug function ax1(ARGS)
function ax1(ARGS)
  tstart=time()

  REMAP={}
  global FS, NFFT, NW, K, PVAL

  local FILELEN_TIC
  local NCHANNELS

  if (length(ARGS)<6)
    require(ARGS[1])
    FILENAME=ARGS[2]
    SUFFIX=ARGS[3]
  else
    FS=ARGS[1]
    NFFT=ARGS[2]
    NW=ARGS[3]
    K=ARGS[4]
    PVAL=ARGS[5]
    FILENAME=ARGS[6]
    SUFFIX=ARGS[7]
  end
  if ((length(ARGS)==5) || (length(ARGS)==9))
    global START
    global STOP
    START=ARGS[end-1]
    STOP=ARGS[end]
  end

  if (isa(FS,String))
    FS = float(FS)
  end
  if (isa(NFFT,String))
    NFFT = int(NFFT)
  end
  if (isa(NW,String))
    NW = int(NW)
  end
  if (isa(K,String))
    K = int(K)
  end
  if (isa(PVAL,String))
    PVAL = float(PVAL)
  end
  if ((length(ARGS)==5) || (length(ARGS)==9))
    if (isa(START,String))
      START = float(START)
    end
    if (isa(STOP,String))
      STOP = float(STOP)
    end
  end

  tmp=splitdir(functionlocs(do_it)[1][1])
  VERSIONAX=readall(`git --git-dir $(tmp[1])/.git log -1 --pretty=format:"%ci %H"`)
  TIMESTAMP=strftime("%Y%m%dT%H%M%S",time())

  SUBSAMPLE=1
  NWORKERS=nworkers()

  FS=int(FS/SUBSAMPLE)

  tmp=log2(NFFT)
  if(abs(tmp-round(tmp))>eps(tmp))
    warn("ax1 will be faster if NFFT is a power of 2")
  end
  tmp=NFFT/2
  if(abs(tmp-round(tmp))>eps(tmp))
    error("NFFT must be even")
  end
  NWINDOWS_PER_WORKER=int(6*256*1000/NFFT)  # NFFT/2 ticks

  tapers = float32(dpss(NFFT,NW))

  f=[0:(NFFT>>1)]*FS/NFFT
  df=f[2]-f[1]

  sig=invlogcdf(FDist(2,2*K-2),log(1-PVAL/NFFT))

  FILEPATH,tmp = splitdir(FILENAME)
  BASEIN,FILETYPE = splitext(tmp)
  FILENAME = joinpath(FILEPATH,BASEIN)
  if isfile(FILENAME*FILETYPE)
    FILEPATH=""
    FILENAMES=[]
    try
      FILELEN_TIC, NCHANNELS=wavread(FILENAME*FILETYPE; format="size")
    catch
      print("can't open file '$FILENAME$FILETYPE'")
      return
    end
    tmp=wavread(FILENAME*FILETYPE; subrange=1)
    if tmp[2]!=FS
      warn("sampling rates in argument list ($FS) and file ($tmp[2]) do not match;  continuing with $tmp[2]")
      FS=tmp[2]
    end
    REMAP=1.0:float(NCHANNELS)
  else
    FILETYPE=".ch"
    tmp = split(readall(`ls $FILEPATH`),"\n")
    tmp2 = map((x) -> ismatch(Regex("$BASEIN.ch[0-9]"),x), tmp)
    FILENAMES=tmp[tmp2]
    NCHANNELS=length(FILENAMES)
    if length(NCHANNELS)==0
      print(["can't find any .ch files with basename '$FILENAME'"])
      return
    end
    for i = 1:NCHANNELS
      filei = string(FILEPATH,"/",FILENAMES[i])
      local fid
      try
        fid = open(filei,"r")
      catch
        print("can't open file '$filei'")
        return
      end
      seekend(fid)
      FILELEN_TIC=int(position(fid)/4)
      close(fid)
      push!(REMAP,float(string(FILENAMES[i][end])))
    end
  end
  FILE_LEN=FILELEN_TIC/FS

  if !isdefined(Main,:START)
    START_TIC=0;
    STOP_TIC=FILELEN_TIC;
  else
    START_TIC=int(START*FS/(NFFT>>1))*(NFFT>>1);
    STOP_TIC=int(STOP*FS);
  end
  tmp=int((STOP_TIC-START_TIC)/(NFFT>>1)-1)
  @printf("Processing %d channels x %.1f min = %d windows = %.1f chunks of data in %s%s\n",
      NCHANNELS, (STOP_TIC-START_TIC)/FS/60, tmp, tmp/NWINDOWS_PER_WORKER, FILENAME, FILETYPE)

  fid = h5open("$FILENAME-$SUFFIX.ax", "w")
  did=d_create(fid,"/hotPixels",datatype(Float32),((1,4),(-1,4)),"chunk",(1024,4))
  attrs(did)["VERSION"] = VERSIONAX
  attrs(did)["TIMESTAMP"] = TIMESTAMP
  attrs(did)["FS"] = FS
  attrs(did)["NFFT"] = NFFT
  attrs(did)["NW"] = NW
  attrs(did)["K"] = K
  attrs(did)["PVAL"] = PVAL
  attrs(did)["FILENAME"] = FILENAME
  attrs(did)["SUFFIX"] = SUFFIX
  if isdefined(Main,:START);  attrs(did)["START"] = START;  end
  if isdefined(Main,:STOP) ;  attrs(did)["STOP"] = STOP  ;  end
  attrs(did)["NWORKERS"] = NWORKERS
  attrs(did)["NWINDOWS_PER_WORKER"] = NWINDOWS_PER_WORKER

  T = [START_TIC : (NFFT>>1)*NWINDOWS_PER_WORKER : STOP_TIC];
  idx = pmap(do_it, [(FILEPATH, FILENAMES, FILENAME, FILETYPE, FILELEN_TIC, NCHANNELS, t,
     NWINDOWS_PER_WORKER, NFFT, PVAL, FS, tapers, sig) for t in T])

  t=1
  k=0
  tloop=time()
  for i in idx
    if (time()-tloop)>10
      @printf("%d sec processed;  %d%% done\n",
          int((T[t]-START_TIC)/FS), int(100*(T[t]-START_TIC)/(STOP_TIC-START_TIC)))
      tloop=time()
    end

#@bp
    sz=size(i,1)
    set_dims!(did,(k+sz,4))
    did[k+(1:sz),1:3]=i[:,1:3]
    did[k+(1:sz),4]=REMAP[i[:,4]]
    k+=sz
    t+=1
  end

  close(fid)

  tstop = time() - tstart
  @printf("Run time was %f minutes.\n",tstop/60)
end

#@iprofile clear
ax1(ARGS)
#@iprofile report
