immutable AxFile
  version::Uint8
  subsample::Uint8
  fs::Uint32
  nfft::Uint32
  nw::Uint16
  k::Uint16
  pval::Float64
  df::Float64
  data::Array{Float64,2}
end

function ==(x::AxFile, y::AxFile)
  if (x.version != y.version)  return false;  end
  if (x.subsample != y.subsample) return false; end
  if (x.fs != y.fs) return false; end
  if (x.nfft != y.nfft) return false; end
  if (x.nw != y.nw) return false; end
  if (x.k != y.k) return false; end
  if (x.pval != y.pval) return false; end
  if (x.df != y.df) return false; end
  if (x.data != y.data) return false; end
  return true
end

function axread(filename)
  fid=open(filename,"r")
  
  VERSION=read(fid,Uint8)
  SUBSAMPLE=read(fid,Uint8)
  read(fid,Uint8)
  FS=read(fid,Uint32)
  NFFT=read(fid,Uint32)
  NW=read(fid,Uint16)
  K=read(fid,Uint16)
  PVAL=read(fid,Float64)
  DF=read(fid,Float64)

  first=position(fid)
  seekend(fid)
  skip(fid,-1)
  last=position(fid)
  if read(fid,Uint8)!='Z'  error("bad .ax file");  end
  seek(fid,first)
  data=transpose(read(fid, Float64, (4,int((last-first)/8/4))))

  close(fid)
  
  AxFile(VERSION,SUBSAMPLE,FS,NFFT,NW,K,PVAL,DF,data)
end
