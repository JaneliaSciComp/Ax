#!/home/arthurb/src/julia/julia

# julia -p NWORKERS ax1.jl params_file FILEIN FILEOUT
# julia -p NWORKERS ax1.jl params_file FILEIN FILEOUT START STOP
# julia -p NWORKERS ax1.jl FS NFFT NW K PVAL FILEIN FILEOUT
# julia -p NWORKERS ax1.jl FS NFFT NW K PVAL FILEIN FILEOUT START STOP
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
# FILEIN: the path and base filename of a single .wav file containing all channels, or
#     or .ch[0-9] files containing arrays of float32
# FILEOUT: an integer to append to FILEIN to differentiate parameter sets used
# START,STOP: optional time range, in seconds
#
# output is a binary file with a time x frequency x amplitude x channel
#     array of hot pixels
#
# julia -p 12 ax1.jl 'ultrasonic_params' 'urine' '1'
# julia -p 12 ax1.jl 200e3 0.001 15 29 0.01 'urine' '1'
# julia -p 12 ax1.jl 450450 0.001 15 29 0.01 0 30 'groundtruth' '1'

#using MAT
using DSP
using WAV
#using Debug
#using Profile
using Distributions

#using ClusterManagers
#addprocs(10, cman=SGEManager())
#require("/home/arthurb/projects/egnor/ax/ax1b.jl")  # full path on cmd line too

require("ax1b.jl")

#@debug
#@bp
function main(ARGS)
  tstart=time()

  REMAP={}
  local FILELEN_TIC
  local NCHANNELS
  local t_now_sec
  local t_offset_tic

  if (length(ARGS)<6)
    require(ARGS[1])
    FILEIN=ARGS[2]
    FILEOUT=ARGS[3]
  else
    FS=ARGS[1]
    NFFT=ARGS[2]
    NW=ARGS[3]
    K=ARGS[4]
    PVAL=ARGS[5]
    FILEIN=ARGS[6]
    FILEOUT=ARGS[7]
  end
  if ((length(ARGS)==5) || (length(ARGS)==9))
    global START
    global STOP
    START=ARGS[end-1]
    STOP=ARGS[end]
  end

  if (isa(FS,String))
    FS = int(FS)
  end
  if (isa(NFFT,String))
    NFFT = float(NFFT)
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

  VERSION=1

  SUBSAMPLE=1
  NWORKERS=nworkers()

  FS=int(FS/SUBSAMPLE)

  NFFT=nextpow2(int(NFFT*FS))  # convert to ticks

  NWINDOWS_PER_WORKER=int(6*256*1000/NFFT)  # NFFT/2 ticks

  tapers = float32(dpss(NFFT,NW))

  f=[0:(NFFT>>1)]*FS/NFFT
  df=f[2]-f[1]

  sig=invlogcdf(FDist(2,2*K-2),log(1-PVAL/NFFT))

  if stat("$FILEIN.wav").inode!=0   # is there a better way?
    FILETYPE="wav"
    FILEPATH=""
    FILEINs=[]
    try
      FILELEN_TIC, NCHANNELS=wavread("$FILEIN.wav"; format="size")
    catch
      print("can't open file '$filei'")
      return
    end
    REMAP=1.0:float(NCHANNELS)
  else
    FILETYPE="ch"
    tmp = split(FILEIN,"/")
    BASEIN = tmp[end]
    FILEPATH = join(tmp[1:end-1],"/")
    tmp = readall(`ls $FILEPATH`)
    tmp = split(tmp,"\n")
    tmp2 = map((x) -> ismatch(Regex("$BASEIN.ch[0-9]"),x), tmp)
    FILEINs=tmp[tmp2]
    NCHANNELS=length(FILEINs)
    if length(NCHANNELS)==0
      print(["can't find any .wav or .ch files with basename '$FILEIN'"])
      return
    end
    for i = 1:NCHANNELS
      filei = string(FILEPATH,"/",FILEINs[i])
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
      push!(REMAP,float(string(FILEINs[i][end])))
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
  @printf("Processing %d channels x %.1f min = %d windows = %.1f chunks of data in %s.%s\n",
      NCHANNELS, (STOP_TIC-START_TIC)/FS/60, tmp, tmp/NWINDOWS_PER_WORKER, FILEIN, FILETYPE)

  fid_out=open("$FILEIN-$FILEOUT.ax","w")
  write(fid_out,uint8(VERSION))
  write(fid_out,uint8(SUBSAMPLE))
  write(fid_out,uint8(0))
  write(fid_out,uint32(FS))
  write(fid_out,uint32(NFFT))
  write(fid_out,uint16(NW))
  write(fid_out,uint16(K))
  write(fid_out,float64(PVAL))
  write(fid_out,float64(df))

  T = [START_TIC : (NFFT>>1)*NWINDOWS_PER_WORKER : STOP_TIC];
  idx = pmap(do_it, [(FILEPATH, FILEINs, FILEIN, FILETYPE, FILELEN_TIC, NCHANNELS, t,
     NWINDOWS_PER_WORKER, NFFT, NW, K, PVAL, FS, tapers, sig) for t in T])

  t=1;
  tloop=time()
  for i in idx
    if (time()-tloop)>10
      @printf("%d sec processed;  %d%% done\n",
          int((T[t]-START_TIC)/FS), int(100*(T[t]-START_TIC)/(STOP_TIC-START_TIC)))
      tloop=time()
    end

    for j in i
      write(fid_out, j[1:3], REMAP[j[4]])
    end
    t+=1
  end

  write(fid_out,"Z")
  close(fid_out)

  tstop = time() - tstart
  @printf("Run time was %f minutes.\n",tstop/60)
end

#@iprofile clear
main(ARGS)
#@iprofile report
