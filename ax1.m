% function ax1(params_file,FILEIN,FILEOUT)
% function ax1(params_file,FILEIN,FILEOUT,START,STOP)
% function ax1(FS,NFFT,NW,K,PVAL,FILEIN,FILEOUT)
% function ax1(FS,NFFT,NW,K,PVAL,FILEIN,FILEOUT,START,STOP)
%
% analyze a set of time series with multi-taper spectral analysis and
% create a sparse matrix of just the time-frequency pixels whose F-test
% passes PVAL.
%
% typical usage consists of one or more input files being analyzed by one
% or more parameter sets.  for example, four microphone recordings of the
% same vocalizing mouse analyzed with three different NFFTs and the same
% NW, K, and PVAL.  <filename>.ch[1-4] yield <filename>-[1-3].ax
%
% FS: sampling rate in Hertz
% NFFT: FFT window size in seconds, rounds up to the next power of 2 tics
% NW: multi-taper time-bandwidth product
% K: number of tapers
% PVAL: F-test p-val threshold
% FILEIN: the base filename and path of .ch[0-9] files containing arrays of float32
% FILEOUT: an integer to append to FILEIN to differentiate parameter sets used
% START,STOP: optional time range, in seconds
%
% output is a binary file with a time x frequency x amplitude x channel array of hot pixels
%
% ax1('ultrasonic_params','urine','1');
% ax1(200e3,0.001,15,29,0.01,'urine','1');
% ax1(450450,0.001,15,29,0.01,0,30,'groundtruth','1');

function ax1(varargin)

if((nargin~=3)&&(nargin~=5)&&(nargin~=7)&&(nargin~=9))
  error('invalid args');
end

close_it=0;
if((exist('matlabpool')==2) && (matlabpool('size')==0))
  try
    matlabpool open
    close_it=1;
  catch
    disp('WARNING: could not open matlab pool.  proceeding with a single thread.');
  end
end

tstart=tic;

if(nargin<6)
  %run(varargin{1});
  fid=fopen(varargin{1},'r');
  eval(fread(fid,'*char')');
  fclose(fid);
  FILEIN=varargin{2};
  FILEOUT=varargin{3};
else
  FS=varargin{1};
  NFFT=varargin{2};
  NW=varargin{3};
  K=varargin{4};
  PVAL=varargin{5};
  FILEIN=varargin{6};
  FILEOUT=varargin{7};
end
if((nargin==5)||(nargin==9))
  START=varargin{end-1};
  STOP=varargin{end};
end

if(ischar(FS))        FS=str2num(FS);              end
if(ischar(NFFT))      NFFT=str2num(NFFT);          end
if(ischar(NW))        NW=str2num(NW);              end
if(ischar(K))         K=str2num(K);                end
if(ischar(PVAL))      PVAL=str2num(PVAL);          end
if((nargin==5)||(nargin==9))
  if(ischar(START))   START=str2num(START);        end
  if(ischar(STOP))    STOP=str2num(STOP);          end
end

if(length(NFFT)>1)
  error('multiple NFFTs not supported when calling ax1() from the matlab command line');
end

VERSION=1;

SUBSAMPLE=1;
NWORKERS=0;
if(exist('matlabpool')==2)
  NWORKERS=matlabpool('size');
end
if(NWORKERS==0)  NWORKERS=1;  end

FS=FS/SUBSAMPLE;

NFFT=2^nextpow2(NFFT*FS);  % convert to ticks
CHUNK=round(12*256*1000/NFFT);  % NFFT/2 ticks

FIRST_MT=nan;
LAST_MT=nan;
FRACTION_MT=nan;

[tapers,eigs]=dpss(NFFT,NW,K);
%tapers = tapers*sqrt(FS);

MT_PARAMS=[];
MT_PARAMS.NW=NW;
MT_PARAMS.K=K;
MT_PARAMS.NFFT=NFFT;
MT_PARAMS.tapers=single(tapers);
MT_PARAMS.Fs=FS;
MT_PARAMS.pad=0;
MT_PARAMS.fpass=[0 FS/2];

f=(0:(NFFT/2))*FS/NFFT;
df=f(2)-f(1);

[FILEPATH n e]=fileparts(FILEIN);
DIR_OUT=fullfile(FILEPATH);
FILEINs=dir([FILEIN '.ch*']);
if(length(FILEINs)==0)
  error(['can''t find file ''' FILEIN '.ch*''']);
end
NCHANNELS=length(FILEINs);

for i=1:length(FILEINs)
  fid=fopen(fullfile(FILEPATH,FILEINs(i).name),'r');
  if(fid==-1)
    error(['can''t open file ''' fullfile(FILEPATH,FILEINs(i).name) '''']);
  end
  fseek(fid,0,1);
  FILE_LEN=ftell(fid)/4/FS;
  if(~exist('START','var'))
    tmp=round(FILE_LEN*FS/(NFFT/2)-1);
    disp(['Processing ' num2str(FILE_LEN/60,3) ' min = ' num2str(tmp) ' windows = ' num2str(tmp/CHUNK,3) ' chunks of data in ' FILEINs(i).name]);
    t_now_offset = 0;
    t_now_sec=0;
  else
    tmp=round((STOP-START)*FS/(NFFT/2)-1);
    disp(['Processing ' num2str((STOP-START)/60,3) ' min = ' num2str(tmp) ' windows = ' num2str(tmp/CHUNK,3) ' chunks of data in ' FILEINs(i).name]);
    t_now_offset = round(START*FS);
    t_now_sec=START;
  end
  REMAP(i)=str2num(FILEINs(i).name(end));
  fclose(fid);
end

fid_out=fopen([FILEIN '-' FILEOUT '.ax'],'w');
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

  idx=cell(1,NWORKERS);
  parfor i=1:NWORKERS
%   for i=1:NWORKERS

    NSAMPLES = NFFT/2*(CHUNK+1);
    dd = zeros(NCHANNELS, NSAMPLES, 'single');
    for j=1:NCHANNELS
      fid = fopen(fullfile(FILEPATH,FILEINs(j).name),'r');
      fseek(fid,((t_now+(i-1)*CHUNK)*NFFT/2+t_now_offset)*4,-1);
      [tmp count] = fread(fid, NSAMPLES, 'float32', 4*(SUBSAMPLE-1));
      if(count<NSAMPLES)
        tmp=[tmp; zeros(NSAMPLES-count, 1, 'single')];
      end
      dd(j,:) = tmp;
      fclose(fid);
    end

    for j=1:CHUNK
      ddd=dd(:,(1:NFFT)+NFFT/2*(j-1));
      [F,p,f,sig,sd] = ftestc(ddd',MT_PARAMS,PVAL/NFFT,'n');
      for l=1:NCHANNELS
        tmp=1+find(F(2:end,l)'>sig);
        for m=1:length(tmp)
          [freq,amp]=brown_puckette(ddd(l,:),f,tmp(m),FS);
          idx{i}{end+1} = [j+(i-1)*CHUNK, freq, amp, l];
        end
      end
    end
  end
  for i=idx
    for j=i{1}
      fwrite(fid_out,[t_now+j{1}(1) j{1}(2:3) REMAP(j{1}(4))],'double');  % blegh
    end
  end

  t_now_sec=t_now_sec+NFFT/2/FS*NWORKERS*CHUNK;
  t_now=t_now+NWORKERS*CHUNK;
end

fwrite(fid_out,'Z','uchar');
fclose(fid_out);

tstop=toc(tstart);
disp(['Run time was ' num2str(tstop/60,3) ' minutes.']);

if((exist('matlabpool')==2) && (matlabpool('size')>0) && close_it)
  try
    matlabpool close
  catch
    disp('WARNING: could not close matlab pool.  exiting anyway.');
  end
end


% from Chronux
function J=mtfftc(data,tapers,nfft,Fs)
if nargin < 4; error('Need all input arguments'); end;
[NC,C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data=data(:,ones(1,K),:); % add taper indices to data
data_proj=data.*tapers; % product of data with tapers
J=fft(data_proj,nfft);   % fft of projected data


% from Chronux
function [Fval,A,f,sig,sd] = ftestc(data,params,p,plt)
if nargin < 1; error('Need data'); end;
if nargin < 2 || isempty(params); params=[]; end;

%[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
tapers=params.tapers;
pad=params.pad;
Fs=params.Fs;
fpass=params.fpass;

% data=change_row_to_column(data);
[N,C]=size(data);
if nargin<3 || isempty(p);p=0.05/N;end;
if nargin<4 || isempty(plt); plt='n';end;

% tapers=dpsschk(tapers,N,Fs); % calculate the tapers
[N,K]=size(tapers);
nfft=max(2^(nextpow2(N)+pad),N);% number of points in fft

%[f,findx]=getfgrid(Fs,nfft,fpass);% frequency grid to be returned
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
findx=find(f>=fpass(1) & f<=fpass(end));
f=f(findx);

Kodd=1:2:K;
Keven=2:2:K;
J=mtfftc(data,tapers,nfft,Fs);% tapered fft of data - f x K x C
Jp=J(findx,Kodd,:); % drop the even ffts and restrict fft to specified frequency grid - f x K x C
H0 = squeeze(sum(tapers(:,Kodd),1)); % calculate sum of tapers for even prolates - K x C 
Nf=length(findx);% number of frequencies
H0sq = sum(H0.*H0)*ones(Nf,C);
H0=H0(ones(1,Nf),:,ones(1,C));
JpH0=squeeze(sum(Jp.*H0,2));% sum of the product of Jp and H0 across taper indices - f x C
A=JpH0./H0sq; % amplitudes for all frequencies and channels
Kp=size(Jp,2); % number of even prolates
Ap=permute(A,[1 3 2]); % permute indices to match those of H0
Ap=Ap(:,ones(1,Kp),:); % add the taper index to C
Jhat=Ap.*H0; % fitted value for the fft

num=(K-1).*(abs(A).^2).*squeeze(H0sq);%numerator for F-statistic
%den=squeeze(sum(abs(Jp-Jhat).^2,2)+sum(abs(J(findx,Keven,:)).^2,2));% denominator for F-statistic
den1=Jp-Jhat;
den1=real(den1).*real(den1)+imag(den1).*imag(den1);
den1=squeeze(sum(den1,2));
den2=J(findx,Keven,:);
den2=real(den2).*real(den2)+imag(den2).*imag(den2);
den2=squeeze(sum(den2,2));
den = den1 + den2;
Fval=num./den; % F-statisitic

sig=finv(1-p,2,2*K-2); % F-distribution based 1-p% point
var=den./(K*squeeze(H0sq)); % variance of amplitude
sd=sqrt(var);% standard deviation of amplitude

A=A*Fs;


% from charpentier (1986) and brown and puckette (1993; JASA)
function [freq,amp]=brown_puckette(x,f,k,fs)

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
