# would be nice to have a subrange like wavread
# should think about missing e.g. .ch1

function chread(filename::String)
  tmp = split(filename,"/")
  BASEIN = tmp[end]
  FILEPATH = join(tmp[1:end-1],"/")
  tmp = readall(`ls $FILEPATH`)
  tmp = split(tmp,"\n")
  tmp2 = map((x) -> ismatch(Regex("$BASEIN.ch[0-9]"),x), tmp)
  FILEINs=tmp[tmp2]

  local data

  for i=1:length(FILEINs)
    fid=open(string(FILEPATH,"/",FILEINs[i]),"r");
    seekend(fid);
    len=position(fid);
    seekstart(fid);
    tmp=read(fid,Float32,int(len/4));
    if(i==1)  data=similar(tmp, (size(tmp,1), length(FILEINs)))  end
    data[:,i]=tmp
    close(fid);
  end

  return data
end
