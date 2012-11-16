% function mtbp(FILEIN,FILEOUT,params_file)
% function mtbp(FILEIN,FILEOUT,params_file,START,STOP)
% function mtbp(FILEIN,FILEOUT,FS,NFFT,NW,K,PVAL)
% function mtbp(FILEIN,FILEOUT,FS,NFFT,NW,K,PVAL,START,STOP)
%
% FS: sampling rate in Hertz
% NFFT: FFT window size in seconds, rounds up to the next power of 2 tics
% NW: multi-taper time-bandwidth product
% K: number of tapers
% PVAL:  F-test p-val threshold
% START,STOP:  optional time range, in seconds
%
% mtbp('urine',1,200e3,0.001,15,29,0.01);
% mtbp('groundtruth',1,450450,0.001,15,29,0.01,0,30);

%function mtbp(FILEIN,FILEOUT,FS,NFFT,NW,K,PVAL,varargin)
function mtbp(FILEIN,FILEOUT,varargin)

if((nargin~=3)&&(nargin~=5)&&(nargin~=7)&&(nargin~=9))
  error('invalid args');
end

tstart=tic;

close_it=0;
if((exist('matlabpool')==2) && (matlabpool('size')==0))
  try
    matlabpool open
    close_it=1;
  catch
    disp('WARNING: could not open matlab pool.  proceeding with a single thread.');
  end
end

if(nargin<6)
  run(varargin{1});
else
  FS=varargin{1};
  NFFT=varargin{2};
  NW=varargin{3};
  K=varargin{4};
  PVAL=varargin{5};
end
if((nargin==5)||(nargin==9))
  START=varargin{end-1};
  STOP=varargin{end};
end

if(ischar(FS))        FS=str2num(FS);              end
if(ischar(NFFT))  NFFT=str2num(NFFT);  end
if(ischar(NW))        NW=str2num(NW);              end
if(ischar(K))         K=str2num(K);                end
if(ischar(PVAL))      PVAL=str2num(PVAL);          end
if((nargin==5)||(nargin==9))
  if(ischar(START))   START=str2num(START);        end
  if(ischar(STOP))    STOP=str2num(STOP);          end
end

if(length(NFFT)>1)
  error('multiple NFFTs not supported when calling mtbp() from the matlab command line');
end

VERSION=1;

SUBSAMPLE=1;
NWORKERS=matlabpool('size');
if(NWORKERS==0)  NWORKERS=1;  end

FS=FS/SUBSAMPLE;

NFFT=2^nextpow2(NFFT*FS);  % convert to ticks
CHUNK=round(256*1000/NFFT);  % NFFT/2 ticks

FIRST_MT=nan;
LAST_MT=nan;
FRACTION_MT=nan;

if ~isdeployed
  %[p n e]=fileparts(which('mtbp'));
  addpath(genpath('~/matlab/chronux'));
end

MT_PARAMS=[];
MT_PARAMS.NW=NW;
MT_PARAMS.K=K;
MT_PARAMS.NFFT=NFFT;
MT_PARAMS.tapers=dpsschk([NW K],NFFT,FS);
MT_PARAMS.Fs=FS;
MT_PARAMS.pad=0;
MT_PARAMS.fpass=[0 FS/2];

f=(0:(NFFT/2))*FS/NFFT;
df=f(2)-f(1);

REMAP=[1:4 6:8];  % blegh

[p n e]=fileparts(FILEIN);
DIR_OUT=fullfile(p);
FILEINs=dir([FILEIN '.ch*']);
[tmp{1:length(FILEINs)}]=deal(FILEINs.name);
tmp=cellfun(@(x) regexp(x,'\.ch[1-4,6-8]'), tmp,'uniformoutput',false);  % blegh
FILEINs=FILEINs(~cellfun(@isempty,tmp));
if(length(FILEINs)==0)
  error(['can''t find file ''' FILEIN '.ch*''']);
end
NCHANNELS=length(FILEINs);

for i=1:length(FILEINs)
  fid(i)=fopen(fullfile(p,FILEINs(i).name),'r');
  if(fid(i)==-1)
    error(['can''t open file ''' fullfile(p,FILEINs(i).name) '''']);
  end
  fseek(fid(i),0,1);
  FILE_LEN=ftell(fid(i))/4/FS;
  if(~exist('START','var'))
    disp(['Processing ' num2str(FILE_LEN/60,3) ' minutes of data in ' FILEINs(i).name]);
    fseek(fid(i),0,-1);
    t_now_sec=0;
  else
    disp(['Processing ' num2str((STOP-START)/60,3) ' minutes of data in ' FILEINs(i).name]);
    fseek(fid(i),round(START*FS)*4,-1);
    t_now_sec=START;
  end
end

dd=zeros(length(fid),NFFT/2*(NWORKERS*CHUNK+1));
for i=1:length(fid)
  dd(i,(end-NFFT/2+1):end)=fread(fid(i),NFFT/2,'float32',4*(SUBSAMPLE-1));
end

fid_out=fopen([FILEIN '-' FILEOUT '.mtbp'],'w');
fwrite(fid_out,uint8([VERSION SUBSAMPLE CHUNK]),'uint8');  % CHUNK not necessary
fwrite(fid_out,uint32([FS NFFT]),'uint32');
fwrite(fid_out,uint16([NW K]),'uint16');
fwrite(fid_out,[PVAL df],'double');

t_now=0;
tic;
while((t_now_sec<FILE_LEN) && (~exist('STOP','var') || (t_now_sec<STOP)))
  if(toc>10)
    tmp=t_now_sec;
    tmp2=0;  if(exist('START','var'))  tmp=tmp-START;  tmp2=START;  end
    if(exist('STOP','var'))  tmp=tmp/(STOP-tmp2);  else  tmp=tmp/(FILE_LEN-tmp2);  end
    disp([num2str(round(t_now_sec-tmp2)) ' sec processed;  ' num2str(round(100*tmp)) '% done']);
    tic;
  end

  dd(:,1:(NFFT/2))=dd(:,(end-NFFT/2+1):end);
  for i=1:length(fid)
    [tmp count]=fread(fid(i),NFFT/2*NWORKERS*CHUNK,'float32',4*(SUBSAMPLE-1));
    if(count<NFFT/2*NWORKERS*CHUNK)
      tmp=[tmp; zeros(NFFT/2*NWORKERS*CHUNK-count,1)];
    end
    dd(i,(NFFT/2+1):end)=tmp;
  end

  idx=cell(NCHANNELS,CHUNK,NWORKERS);
  parfor i=1:NWORKERS
    for j=1:CHUNK
      [F,p,f,sig,sd] = ftestc(dd(:,(1:NFFT)+NFFT/2*(j+(i-1)*CHUNK-1))',MT_PARAMS,PVAL/NFFT,'n');
      for l=1:NCHANNELS
        tmp=1+find(F(2:end,l)'>sig);
        tmp2=[];
        for m=1:length(tmp)
          [tmp2(m,1) tmp2(m,2)]=brown2_puckette(dd(l,(1:NFFT)+NFFT/2*(j+(i-1)*CHUNK-1)),f,tmp(m),FS);
        end
        idx{l,j,i}=tmp2;
      end
    end
  end
  idx=reshape(idx,NCHANNELS,NWORKERS*CHUNK);
  [sub1,sub2]=ind2sub(size(idx),find(~cellfun(@isempty,idx)));
  if(length(sub1)>0)
    for i=1:length(sub1)
      tmp=idx{sub1(i),sub2(i)};
      for j=1:size(tmp,1)
        fwrite(fid_out,[t_now+sub2(i) tmp(j,1) tmp(j,2) REMAP(sub1(i))],'double');  % blegh
      end
    end
  else
  end

  t_now_sec=t_now_sec+NFFT/2/FS*NWORKERS*CHUNK;
  t_now=t_now+NWORKERS*CHUNK;
end

fwrite(fid_out,'Z','uchar');
fclose(fid_out);
for i=1:length(FILEINs)
  fclose(fid(i));
end

tstop=toc(tstart);
disp(['Run time was ' num2str(tstop/60,3) ' minutes.']);

if((exist('matlabpool')==2) && (matlabpool('size')>0) && close_it)
  try
    matlabpool close
  catch
    disp('WARNING: could not close matlab pool.  exiting anyway.');
  end
end


function [freq,amp]=brown2_puckette(x,f,k,fs)

nfft=length(x);
X=fft(x);
Xh0=0.5*(X(k)-0.5*X(k+1)-0.5*X(k-1));
Xh1=0.5*exp(sqrt(-1)*2*pi*(k-1)/nfft)*...
   (X(k) - 0.5*exp(sqrt(-1)*2*pi/nfft)*X(k+1)...
         - 0.5*exp(-sqrt(-1)*2*pi/nfft)*X(k-1));
phi0=atan2(imag(Xh0),real(Xh0));
phi1=atan2(imag(Xh1),real(Xh1));
if((phi1-phi0)<0)  phi1=phi1+2*pi;  end
freq=(phi1-phi0)*fs/(2*pi);

period = fs/freq;
last = floor(period * floor(length(x)/period));
real_part = mean(x(1:last) .* cos([1:last]*(2*pi/period)));
imag_part = mean(x(1:last) .* sin([1:last]*(2*pi/period)));
amp = 2*abs(real_part + i*imag_part);
