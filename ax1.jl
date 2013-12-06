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
# FILEIN: the base filename and path of [0-9].wav files with a single channel each,
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

require("ax1b.jl")
using MAT
using DSP
using WAV
#using Debug
#using Profile

function main(ARGS)
  tstart=time()

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

  NWINDOWS_PER_WORKER=int(12*256*1000/NFFT)  # NFFT/2 ticks

  FIRST_MT=NaN
  LAST_MT=NaN
  FRACTION_MT=NaN

  f=[0:(NFFT>>1)]*FS/NFFT
  df=f[2]-f[1]

  tmp = split(FILEIN,"/")
  BASEIN = tmp[end]
  DIROUT = join(tmp[1:end-1],"/")
  tmp = readall(`ls $DIROUT`)
  tmp = split(tmp,"\n")
  tmp2 = map((x) -> ismatch(Regex("$BASEIN.ch[0-9]"),x), tmp)
  FILE_TYPE=1
  if !any(tmp2)
    tmp2 = map((x) -> ismatch(Regex("$BASEIN[0-9].wav"),x), tmp)
    FILE_TYPE=2
    if !any(tmp2)
      print(["can't find any .wave or .ch files with basenmae '$FILEIN'"])
      return
    end
  end
  FILEINs=tmp[tmp2]
  NCHANNELS=length(FILEINs)

  REMAP={}
  local FILE_LEN
  local t_now_sec
  local t_offset_tic
  for i = 1:NCHANNELS
    filei = string(DIROUT,"/",FILEINs[i])
    if FILE_TYPE==1
      local fid
      try
        fid = open(filei,"r")
      catch
        print("can't open file '$filei'")
        return
      end
      seekend(fid)
      FILE_LEN=position(fid)/4/FS
      close(fid)
      push!(REMAP,float(string(FILEINs[i][end])))
    end
    if FILE_TYPE==2
      try
        FILE_LEN, foo = wavread(filei,format="size")
      catch
        print("can't open file '$filei'")
        return
      end
      FILE_LEN /= FS
      push!(REMAP,float(string(FILEINs[i][end-4])))
    end
    if !isdefined(Main,:START)
      tmp=int(FILE_LEN*FS/(NFFT>>1)-1)
      @printf("Processing %f min = %d windows = %f chunks of data in %s\n",
          FILE_LEN/60, tmp, tmp/NWINDOWS_PER_WORKER, FILEINs[i])
      t_offset_tic=0;
      t_now_sec=0;
    else
      tmp=int((STOP-START)*FS/(NFFT>>1)-1)
      @printf("Processing %f min = %d windows = %f chunks of data in %s\n",
          (STOP-START)/60, tmp, tmp/NWINDOWS_PER_WORKER, FILEINs[i])
      t_offset_tic=int(START*FS);
      t_now_sec=START;
    end
  end

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

  tapers = float32(dpss(NFFT,NW))

  t_now=0
  tloop=time()

  while ((t_now_sec<FILE_LEN) && (!isdefined(Main,:STOP) || (t_now_sec<STOP)))
    if (time()-tloop)>10
      tmp=t_now_sec
      tmp2=0
      if isdefined(Main,:START)
        tmp=tmp-START
        tmp2=START
      end
      if isdefined(Main,:STOP)
        tmp=tmp/(STOP-tmp2)
      else
        tmp=tmp/(FILE_LEN-tmp2)
      end
      @printf("%d sec processed;  %d%% done\n",int(round(t_now_sec-tmp2)),int(round(100*tmp)))
      tloop=time()
    end

    idx = pmap(do_it,
       [(DIROUT, FILEINs, t_now, NW, K, PVAL, FS, NFFT, NWINDOWS_PER_WORKER, tapers, x-1, t_offset_tic, FILE_TYPE, int(round(FILE_LEN*FS))) for x in 1:NWORKERS])

    for i in idx
      for j in i
        write(fid_out, t_now+j[1], j[2:3], REMAP[j[4]])
      end
    end

    t_now_sec = t_now_sec+float(NFFT>>1)/FS*NWORKERS*NWINDOWS_PER_WORKER
    t_now = t_now+NWORKERS*NWINDOWS_PER_WORKER
  end

  write(fid_out,"Z")
  close(fid_out)

  tstop = time() - tstart
  @printf("Run time was %f minutes.\n",tstop/60)
end

#addprocs(12, cman=SGEManager())

#@iprofile clear
main(ARGS)
#@iprofile report
