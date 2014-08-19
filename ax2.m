% function directory=ax2(frequency_low, frequency_high, convolution_size, minimum_object_area, ...
%     merge_harmonics, merge_harmonics_overlap, merge_harmonics_ratio, merge_harmonics_fraction, ...
%     minimum_vocalization_length, ...
%     channels, data_path)
% function directory=ax2(params_file, data_path)
%
% extract contours from a set of spectrograms by first convolving hot pixels
% with a square box, then finding contiguous pixels, and finally discarding those
% with small areas.  optionally merge contours that are harmonically related.
%
% typical usage consists of three spectrograms made with three different FFT
% window sizes.  a single common image is constructed using the smallest temporal
% and frequency resolutions and populating with horizontally or vertically
% elongated pixels.
%
% input arguments:
%   frequency_low, frequency_high are in Hz
%   convolution_size is [frequency time] in Hz and sec
%   minimum_object_area is in Hz-sec
%   set merge_harmonics to 1 to collapse harmonically related syllables, 0 otherwise
%     merge_harmonics_overlap is the fraction in time two segments must overlap
%     merge_harmonics_ratio is the tolerance in frequency ratio two segments must be within
%     merge_harmonics_fraction is the fraction of the overlap that must be within the ratio tolerance
%   minimum_vocalization_length is the minimum vocalization length in sec
%   channels is a vector of which channels to use, or [] to use all of them
%   data_path can be to a folder or to a set of files.  for the latter, omit the .wav suffix
%
% four files are output:
%   voc: an Mx4 array whose columns are the start & stop times (sec), and low & high frequences (Hz)
%   fc: a cell array (vocs) of cell arrays (syls) of Nx3 arrays (time[s], freq[Hz], amplitude)
%       of the hot pixel in each time slice with the max amplitude
%   fc2: a cell array (vocs) of cell arrays (syls) of Nx4 arrays (time[s], freq[Hz], amplitude, channel)
%       of all hot pixels
%   params: a .m file of the parameters used
%
% for ultrasonic:
% ax2(20e3, 120e3, [15 7], 1500, 0, 0.9, 0.1, 0.9, 0, 1, 0,...
%     1:4, '/groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/');
% ax2('ultrasonic_params.m', '/groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/');
%
% for rejection:
% ax2(1e3, 20e3, [15 7], 1000, 0, 0.9, 0.1, 0.9, 0, 3, 0,...
%     1:4, '/groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/');
% ax2('rejection_params.m', '/groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/');

function ax2(varargin)

switch nargin
  case 2
%    run(varargin{1});
    fid = fopen(varargin{1});
    if fid < 0
        error('Could not open the parameters file at %s', varargin{1});
    end
    params_code = fread(fid, '*char')';
    fclose(fid);
    try
%        disp(params_code);
        eval(params_code);
    catch ME
        error('Could not load the parameters from %s (%s)', varargin{1}, ME.message);
    end
    data_path=varargin{2};
  case 11
    frequency_low=varargin{1};
    frequency_high=varargin{2};
    convolution_size=varargin{3};
    minimum_object_area=varargin{4};
    merge_harmonics=varargin{5};
    merge_harmonics_overlap=varargin{6};
    merge_harmonics_ratio=varargin{7};
    merge_harmonics_fraction=varargin{8};
    minimum_vocalization_length=varargin{9};
    channels=varargin{10};
    data_path=varargin{11};
  otherwise
    error([num2str(nargin) ' args input, either 2 or 11 expected']);
end

if(ischar(frequency_low))                 frequency_low=str2num(frequency_low);                               end
if(ischar(frequency_high))                frequency_high=str2num(frequency_high);                             end
if(ischar(convolution_size))              convolution_size=str2num(convolution_size);                         end
if(ischar(minimum_object_area))           minimum_object_area=str2num(minimum_object_area);                   end
if(ischar(merge_harmonics))               merge_harmonics=str2num(merge_harmonics);                           end
if(ischar(merge_harmonics_overlap))       merge_harmonics_overlap=str2num(merge_harmonics_overlap);           end
if(ischar(merge_harmonics_ratio))         merge_harmonics_ratio=str2num(merge_harmonics_ratio);               end
if(ischar(merge_harmonics_fraction))      merge_harmonics_fraction=str2num(merge_harmonics_fraction);         end
if(ischar(minimum_vocalization_length))   minimum_vocalization_length=str2num(minimum_vocalization_length);   end
if(ischar(channels))                      channels=str2num(channels);                                         end

if(isempty(frequency_low) || isempty(frequency_high) || (frequency_low<0) || (frequency_high<0) || (frequency_low>=frequency_high))
  error('frequency_low should be less than frequency_high and both should be non-negative real numbers');
end

if(isempty(convolution_size) || (length(convolution_size)~=2))
  error('convolution_size should be a 2-vector');
end

if(isempty(minimum_object_area) || (minimum_object_area<0))
  error('minimum_object_area should be a non-negative integer');
end

if (isempty(merge_harmonics) || ((merge_harmonics~=0) && (merge_harmonics~=1)))
  error('merge_harmonics must be 0 or 1');
end

if (isempty(merge_harmonics_overlap) || ((merge_harmonics_overlap<0) || (merge_harmonics_overlap>1)))
  error('merge_harmonics_overlap must be between 0 and 1');
end

if (isempty(merge_harmonics_ratio) || (merge_harmonics_ratio<0))
  error('merge_harmonics_ratio must be a non-negative real number');
end

if (isempty(merge_harmonics_fraction) || (merge_harmonics_fraction<0) || (merge_harmonics_fraction>1))
  error('merge_harmonics_fraction must be between 0 and 1');
end

if (isempty(minimum_vocalization_length) || (minimum_vocalization_length<0))
  warndlg('minimum_vocalization_length must be a non-negative real number');
end

tmp=dir(fullfile(data_path,'-*.ax'));
if(~isempty(tmp))
  datafiles=cell(1,length(tmp));
  for i=1:length(tmp)
    [~,datafiles{i},~]=fileparts(tmp(i).name(1:end-5));
  end
  datafiles=unique(datafiles);
  close_it=0;
  if(~isdeployed && length(datafiles)>1)
    if(exist('parpool')==2 && (length(gcp)==0))
      try
        parpool;
        close_it=1;
      catch
        disp('WARNING: could not open parallel pool of workers.  proceeding with a single thread.');
      end
    end
  end
  parfor i=1:length(datafiles)
%   for i=1:length(datafiles)
    ax2_guts(i, frequency_low, frequency_high, convolution_size, minimum_object_area, ...
        merge_harmonics, merge_harmonics_overlap, merge_harmonics_ratio, merge_harmonics_fraction, ...
        minimum_vocalization_length, channels, fullfile(data_path, datafiles{i}));
  end
  if((exist('parpool')==2) && (length(gcp)==0) && close_it)
    try
      delete(gcp('nocreate'));
    catch
      disp('WARNING: could not close parallel pool of workers.  exiting anyway.');
    end
  end
else
  tmp=dir([data_path '-*.ax']);
  if(~isempty(tmp))
    ax2_guts(0, frequency_low, frequency_high, convolution_size, minimum_object_area, ...
        merge_harmonics, merge_harmonics_overlap, merge_harmonics_ratio, merge_harmonics_fraction, ...
        minimum_vocalization_length, channels, data_path);
  else
    error(['can''t find ' data_path]);
  end
end


function ax2_guts(num, FREQUENCY_LOW, FREQUENCY_HIGH, CONVOLUTION_SIZE, MINIMUM_OBJECT_AREA, ...
    MERGE_HARMONICS, MERGE_HARMONICS_OVERLAP, MERGE_HARMONICS_RATIO, MERGE_HARMONICS_FRACTION, ...
    MINIMUM_VOCALIZATION_LENGTH, CHANNELS, filename)

SAVE_WAV=0;
SAVE_PNG=0;

if(SAVE_WAV || SAVE_PNG)
  figure;
  get(gcf,'position');
  set(gcf,'position',[ans(1) ans(2) 4*ans(3) ans(4)]);
  subplot('position',[0.05 0.1 0.9 0.8]);
  set(gca,'color',[0 0 0]);
end

%load header
disp([num2str(num) ': ' 'processing file ' filename]);
tmp=dir([filename,'-*.ax']);
filenames={tmp.name};
for i=1:length(filenames)
  disp([num2str(num) ': loading ' filenames{i}]);
  try
    hdf5=1;
    tmp2=h5info(fullfile(fileparts(filename),filenames{i}),'/hotPixels');
    for j=1:length(tmp2.Attributes)
      data(i).(tmp2.Attributes(j).Name)=tmp2.Attributes(j).Value;
    end
    flen(i)=tmp2.Dataspace.Size(1);
    fptr(i)=1;
  catch
    hdf5=0;
    fid(i)=fopen(fullfile(fileparts(filename),filenames{i}));
    data(i).VERSION_FILE_FORMAT=fread(fid(i),1,'uint8');
    data(i).SUBSAMPLE=fread(fid(i),1,'uint8');
    data(i).CHUNK=fread(fid(i),1,'uint8');
    data(i).FS=fread(fid(i),1,'uint32');
    data(i).NFFT=fread(fid(i),1,'uint32');
    data(i).NW=fread(fid(i),1,'uint16');
    data(i).K=fread(fid(i),1,'uint16');
    data(i).PVAL=fread(fid(i),1,'double');
    fread(fid(i),1,'double');  % skip deltaFreq
    len=ftell(fid(i));
    fseek(fid(i),-1,'eof');
    if(~strcmp(char(fread(fid(i),1,'uchar')),'Z'))
      error('end of file marker missing.  file corrupt');
    end
    fseek(fid(i),len,'bof');
  end
  if((i>1) && (data(i).FS~=data(1).FS))
    error('sampling frequencies are not the same');
  end
end
[~,idx]=sort([data.NFFT]);
data=data(idx);
filenames=filenames(idx);
if(hdf5)
  flen=flen(idx);
else
  fid=fid(idx);
end

deltaFreq=min([data.FS]./[data.NFFT])/10;
minNFFT=min([data.NFFT]);
maxNFFT=max([data.NFFT]);
FS=data(1).FS;
CHUNK_TIME_SEC=5;  % sec
CHUNK_TIME_WINDOWS=round(CHUNK_TIME_SEC*FS/(maxNFFT/2))*maxNFFT./[data.NFFT];  % in units of windows
skytruth=[];
freq_contours={};
freq_contours2={};
voc_num=1;
hit_num=1;
miss_num=1;
fa_num=1;
CHUNK_FILE=1024;
CONVOLUTION_SIZE_PIX=round(CONVOLUTION_SIZE ./ [deltaFreq minNFFT/FS/2]);
CONVOLUTION_SIZE_PIX=round(2*(0.5+floor(CONVOLUTION_SIZE_PIX/2)));  % ceil to odd
MINIMUM_OBJECT_AREA_PIX=MINIMUM_OBJECT_AREA/deltaFreq/(minNFFT/FS/2);

tic;
eof=false;  count=4*CHUNK_FILE;
chunk_curr=1;
while ~eof
  if(toc>10)  disp([num2str(num) ': ' num2str(chunk_curr*CHUNK_TIME_SEC) ' sec chunk']);  tic;  end;

  eof=(max(count)<4*CHUNK_FILE);

  % read in chunk of data
  for i=1:length(data)
    tmp=[];  data(i).MT_next=[];
    while (hdf5 && fptr(i)<=flen(i)) || ((~hdf5) && (~feof(fid(i))))
      if(hdf5)
        count(i)=min(flen(i)-fptr(i)+1, CHUNK_FILE);
        foo=h5read(fullfile(fileparts(filename),filenames{i}),'/hotPixels',[fptr(i) 1],[count(i) 4]);
        fptr(i)=fptr(i)+count(i);
        count(i)=count(i)*4;
      else
        [foo,count(i)]=fread(fid(i),[4 CHUNK_FILE],'double');
        foo=foo';
      end
      if(isempty(foo))  continue;  end
      tmp=[tmp; foo];
      idx=find(tmp(:,1)>chunk_curr*CHUNK_TIME_WINDOWS(i),1);
      if(((hdf5 && fptr(i)>flen(i)) || ((~hdf5) && feof(fid(i)))) && isempty(idx))
        idx=size(tmp,1)+1;
      end
      if(~isempty(idx))
        idx2=find((tmp(1:(idx-1),2)>=FREQUENCY_LOW) & (tmp(1:(idx-1),2)<=FREQUENCY_HIGH) & ...
            (isempty(CHANNELS) | ismember(tmp(1:(idx-1),4),CHANNELS)));
        data(i).MT_next=tmp(idx2,:);
        data(i).MT_next(:,1)=data(i).MT_next(:,1)-(chunk_curr-1)*CHUNK_TIME_WINDOWS(i);
        if(hdf5)
          fptr(i)=fptr(i)-(size(tmp,1)-idx+1);
        else
          fseek(fid(i),-(size(tmp,1)-idx+1)*4*8,'cof');
        end
        break;
      end
    end
  end

  sizeF=ceil(FREQUENCY_HIGH/deltaFreq)-floor(FREQUENCY_LOW/deltaFreq)+2*floor(maxNFFT/minNFFT/2)+1;
  sizeT=CHUNK_TIME_WINDOWS(1)+2*floor(maxNFFT/minNFFT/2)+1;
  %skip this chunk if no hot pixels
  if(~all(arrayfun(@(x) isempty(x.MT_next), data)))
    im_next=false(sizeF,sizeT);

    %collapse across channels and window sizes
    for k=1:length(data)
      if(isempty(data(k).MT_next))  continue;  end
      tmpT=data(k).NFFT/minNFFT;
      tmpF=maxNFFT/data(k).NFFT;
      for i=(-floor(tmpF/2):floor(tmpF/2))+floor(maxNFFT/minNFFT/2)+1
        for n=(-floor(tmpT/2):floor(tmpT/2))+floor(maxNFFT/minNFFT/2)+1
          im_next(sub2ind([sizeF,sizeT],...
              round(data(k).MT_next(:,2)/deltaFreq)+i-floor(FREQUENCY_LOW/deltaFreq), ...
              tmpT*data(k).MT_next(:,1)+n))=true;
        end
      end
    end

    %convolve
    im_next=[false(sizeF,(CONVOLUTION_SIZE_PIX(2)-1)/2) im_next false(sizeF,(CONVOLUTION_SIZE_PIX(2)-1)/2)];
    im_next=logical(conv2(single(im_next),ones(CONVOLUTION_SIZE_PIX),'same'));

    %segment
    syls_next=bwconncomp(im_next,8);
  else
    syls_next.Connectivity=0;
    syls_next.ImageSize=[sizeF sizeT];
    syls_next.NumObjects=0;
    syls_next.PixelIdxList={};
  end

  % skip to 2nd chunk if currently on first
  if(exist('syls'))

    %unsplit across chunk boundaries
    flag=1;
    while(flag)
      flag=0;
      for i=1:syls.NumObjects
        [ri ci]=ind2sub(syls.ImageSize,syls.PixelIdxList{i});
        if(sum(ci>=(syls.ImageSize(2)-(CONVOLUTION_SIZE_PIX(2)-1)/2))==0) continue;  end
        j=1;
        while j<=syls_next.NumObjects
          [rj cj]=ind2sub(syls_next.ImageSize,syls_next.PixelIdxList{j});
          if(sum(cj<=((CONVOLUTION_SIZE_PIX(2)-1)/2))==0)  j=j+1;  continue;  end
          %if((max(ri) < (min(rj)-1)) || (max(rj) < (min(ri)-1)))  j=j+1;  continue;  end
          %cj=cj+chunk_splits(k+1)-chunk_splits(k);
          cj=cj+CHUNK_TIME_WINDOWS(1);
%           min(min((repmat(ri,1,length(rj))-repmat(rj',length(ri),1)).^2 + ...
%                   (repmat(ci,1,length(cj))-repmat(cj',length(ci),1)).^2));
          smallest=inf;
          for k=1:length(rj)
            smallest=min(smallest, min((ri-rj(k)).^2 + (ci-cj(k)).^2));
          end
          if smallest<=2
            %disp(['unsplitting syllable between chunks #' num2str(k) '-' num2str(k+1)]);
            flag=1;
            syls.PixelIdxList{i}=[syls.PixelIdxList{i}; ...
                syls_next.PixelIdxList{j}+CHUNK_TIME_WINDOWS(1)*syls.ImageSize(1)];
            syls_next.PixelIdxList(j)=[];
            syls_next.NumObjects=syls_next.NumObjects-1;
          else
            j=j+1;
          end
        end
      end
    end

    %throwout out small blobs
    syls2=regionprops(syls,'basic');
    tmp=find([syls2.Area]>MINIMUM_OBJECT_AREA_PIX);
    syls2=syls2(tmp);
    syls.PixelIdxList=syls.PixelIdxList(tmp);
    syls.NumObjects=length(tmp);

    %calculate frequency contours
    freq_contour={};
    freq_contour2={};
    for i=1:length(syls2)
      tmp=[];
      for j=1:length(data)
        if(isempty(data(j).MT))  continue;  end
        tmpT=data(j).NFFT/minNFFT;
  %      idx=find(((tmpT*data(j).MT(:,1)+(CONV_SIZE(2)-1)/2-floor(tmpT/2))>=syls2(i).BoundingBox(1)) & ...
  %               ((tmpT*data(j).MT(:,1)+(CONV_SIZE(2)-1)/2+floor(tmpT/2))<=sum(syls2(i).BoundingBox([1 3]))) & ...
  %                (data(j).MT(:,2)>=(syls2(i).BoundingBox(2)*deltaFreq+F_LOW)) & ...
  %                (data(j).MT(:,2)<=(sum(syls2(i).BoundingBox([2 4]))*deltaFreq+F_LOW)));
        foo=[tmpT*data(j).MT(:,1)+floor(maxNFFT/minNFFT/2)+1+(CONVOLUTION_SIZE_PIX(2)-1)/2 ...
             round(data(j).MT(:,2)/deltaFreq)-floor(FREQUENCY_LOW/deltaFreq)+floor(maxNFFT/minNFFT/2)+1];
        [r c]=ind2sub(syls.ImageSize,syls.PixelIdxList{i});
        idx=ismember(foo,[c r],'rows');
        foo=data(j).MT(idx,1:4);
        foo(:,1)=foo(:,1)+(chunk_curr-2)*CHUNK_TIME_WINDOWS(j);  % +1?
        foo(:,1)=foo(:,1).*data(j).NFFT/2/FS;
        tmp=[tmp; foo];
      end
      tmp=sortrows(tmp);
      freq_contour{i}{1}=zeros(length(unique(tmp(:,1))),3);
      freq_contour2{i}{1}=tmp;
      j=1;  l=1;
      while(j<=size(tmp,1))
        k=j+1;  while((k<=size(tmp,1)) && (tmp(j,1)==tmp(k,1)))  k=k+1;  end
        [~, idx]=max(tmp(j:(k-1),3));
        freq_contour{i}{1}(l,:)=tmp(j+idx-1,1:3);
        j=k;  l=l+1;
      end
    end

    %merge harmonically related syllables
    if(MERGE_HARMONICS~=0)
      syls3=ones(1,length(syls2));
      for i=1:(length(syls2)-1)
        if(isempty(syls.PixelIdxList{i}))  continue;  end
        flag=1;
        while(flag)
          flag=0;
          for j=(i+1):length(syls2)
            if(isempty(syls.PixelIdxList{j}))  continue;  end
            doit=false;  position=nan(1,length(freq_contour{i}));
            for k=1:length(freq_contour{i})
              [c,ii,jj]=intersect(freq_contour{i}{k}(:,1),freq_contour{j}{1}(:,1));
              if((length(c)/size(freq_contour{i}{k},1)<MERGE_HARMONICS_OVERLAP) && ...
                 (length(c)/size(freq_contour{j}{1},1)<MERGE_HARMONICS_OVERLAP))
                continue;
              end
              %doit=doit | ...
              %    sum(sum(abs((freq_contour{i}{k}(ii,2)./freq_contour{j}{1}(jj,2))*...
              %        [1/3 1/2 2/3 3/2 2 3]-1)<MERGE_FREQ_RATIO)...
              %    >(MERGE_FREQ_FRACTION*length(c)));
              sum(abs((freq_contour{i}{k}(ii,2)./freq_contour{j}{1}(jj,2))*...
                  [1/3 1/2 2/3 3/2 2 3]-1)<MERGE_HARMONICS_RATIO)>(MERGE_HARMONICS_FRACTION*length(c));
              doit=doit | sum(ans);
              find(ans,1,'first');
              if(~isempty(ans))
                position(k)=ans<4;
              end
            end
            if(doit)
              flag=1;
              syls.PixelIdxList{i}=[syls.PixelIdxList{i}; syls.PixelIdxList{j}];
              syls.PixelIdxList{j}=[];
              syls2(i).BoundingBox(1)=...
                  min([syls2(i).BoundingBox(1) syls2(j).BoundingBox(1)]);
              syls2(i).BoundingBox(2)=...
                  min([syls2(i).BoundingBox(2) syls2(j).BoundingBox(2)]);
              syls2(i).BoundingBox(3)=...
                  max([sum(syls2(i).BoundingBox([1 3])) sum(syls2(j).BoundingBox([1 3]))])-...
                  syls2(i).BoundingBox(1);
              syls2(i).BoundingBox(4)=...
                  max([sum(syls2(i).BoundingBox([2 4])) sum(syls2(j).BoundingBox([2 4]))])-...
                  syls2(i).BoundingBox(2);
              syls3(i)=syls3(i)+syls3(j);
              syls3(j)=0;
              tmp=find(position,1,'first');  if(isempty(tmp))  tmp=length(freq_contour{i})+1;  end
              freq_contour{i}={freq_contour{i}{1:(tmp-1)} freq_contour{j}{:} freq_contour{i}{tmp:end}};
              freq_contour{j}=[];
              freq_contour2{i}={freq_contour2{i}{1:(tmp-1)} freq_contour2{j}{:} freq_contour2{i}{tmp:end}};
              freq_contour2{j}=[];
            end
          end
        end
      end
      idx=find(syls3>=1);
      syls.NumObjects=length(idx);
      syls.PixelIdxList={syls.PixelIdxList{idx}};
      syls2=regionprops(syls,'basic');
      freq_contour={freq_contour{idx}};
      freq_contour2={freq_contour2{idx}};
    end

    %cull short vocalizations
    if(MINIMUM_VOCALIZATION_LENGTH>0)
      reshape([syls2.BoundingBox],4,length(syls2))';
      idx=find(((ans(:,3)-CONVOLUTION_SIZE_PIX(2)+1)*minNFFT/2/FS)>MINIMUM_VOCALIZATION_LENGTH);
      syls.NumObjects=length(idx);
      syls.PixelIdxList={syls.PixelIdxList{idx}};
      syls2=regionprops(syls,'basic');
      freq_contour={freq_contour{idx}};
      freq_contour2={freq_contour2{idx}};
    end
    tmp=reshape([syls2.BoundingBox],4,length(syls2))';
    %tmp(:,1)=tmp(:,1)-(CONV_SIZE(2)-1)/2+(CONV_SIZE(2)-1)/2;
    tmp(:,1)=tmp(:,1)-floor(maxNFFT/minNFFT/2)-1;
    tmp(:,2)=tmp(:,2)-floor(maxNFFT/minNFFT/2)-1;
    tmp(:,3)=tmp(:,3)-CONVOLUTION_SIZE_PIX(2)+1;
    skytruth=[skytruth; ...
        ([tmp(:,1) tmp(:,1)+tmp(:,3)]+(chunk_curr-2)*CHUNK_TIME_WINDOWS(1))*minNFFT/2/FS ...
        zeros(size(tmp,1),1) ...
        [tmp(:,2) tmp(:,2)+tmp(:,4)].*deltaFreq+FREQUENCY_LOW];
    freq_contours={freq_contours{:} freq_contour{:}};
    freq_contours2={freq_contours2{:} freq_contour2{:}};

    %plot
    if(SAVE_WAV || SAVE_PNG)
      clf;  hold on;

      [r,c]=ind2sub(syls.ImageSize,cat(1,syls.PixelIdxList{:}));
      c=c-(CONVOLUTION_SIZE_PIX(2)-1)/2+(chunk_curr-2)*CHUNK_TIME_WINDOWS(1);
      c=c-floor(maxNFFT/minNFFT/2)-1;
      r=r-floor(maxNFFT/minNFFT/2)-1;
      plot(c.*(minNFFT/2)./FS,r.*deltaFreq+FREQUENCY_LOW,'bo');

      for i=1:length(syls2)
        for j=1:length(freq_contours{end-i+1})
          plot(freq_contours{end-i+1}{j}(:,1),freq_contours{end-i+1}{j}(:,2),'r-');
          plot(freq_contours2{end-i+1}{j}(:,1),freq_contours2{end-i+1}{j}(:,2),'g.');
        end
      end

      plot(skytruth((end-length(syls2)+1):end,[1 2 2 1 1])',...
           skytruth((end-length(syls2)+1):end,[4 4 5 5 4])','y');

      %plot(repmat(chunk_splits*minNFFT/2/FS,2,1),...
      %     repmat([F_LOW; F_HIGH],1,length(chunk_splits)),'c');

      axis tight;
      v=axis;  axis([v(1) v(2) FREQUENCY_LOW FREQUENCY_HIGH]);
      xlabel('time (s)');
      ylabel('frequency (Hz)');

      [b,a]=butter(4,FREQUENCY_LOW/(FS/2),'high');

      % plot hits, misses and false alarms separately
      while voc_num<=min([200 size(skytruth,1)])
        left=skytruth(voc_num,1);
        right=skytruth(voc_num,2);
        ax2_print(voc_num,left,right,'voc',filename,FS,minNFFT,SAVE_WAV,SAVE_PNG,b,a,CHANNELS);
        voc_num=voc_num+1;
      end
    end
  end

  for i=1:length(data)
    clear data(i).MT;  data(i).MT=data(i).MT_next;
  end
%   im=im_next;
  clear syls;  syls=syls_next;
  chunk_curr=chunk_curr+1;
end

%dump files
disp([num2str(num) ': ' num2str(size(skytruth,1)) ' automatically segmented vocalizations']);

tmp=[skytruth(:,1:2) skytruth(:,4:5)];
save([filename '.voc'],'tmp','-ascii');
save([filename '.fc'],'freq_contours');
save([filename '.fc2'],'freq_contours2');

varname=@(x) inputname(1);
%fid=fopen(fullfile(directory,['params' sprintf('%d',CHANNELS) '.m']),'w');
fid=fopen([filename '.params'],'w');
for i=1:length(data)
  jj=fieldnames(data(i));
  for j=1:length(jj)-2
    if(strcmp(jj(j),'NFFT'))
      fprintf(fid,'%s=%g;  %% deltaT=%g secs; deltaF=%g Hz\n',...
          char(jj(j)),data(i).(char(jj(j))), ...
          data(i).(char(jj(j))) / data(i).FS, ...
          data(i).FS / data(i).(char(jj(j))));
    elseif(isnumeric(data(i).(char(jj(j)))))
      fprintf(fid,'%s=%g;\n',char(jj(j)),data(i).(char(jj(j))));
    else
      fprintf(fid,'%s=%s;\n',char(jj(j)),data(i).(char(jj(j))));
    end
  end
  fprintf(fid,'\n');
end

VERSION=get_version();

fprintf(fid,'VERSION=''%s'';\n',VERSION);
fprintf(fid,'TIME_STAMP=''%s'';\n',datestr(now,30));
fprintf(fid,'%s=%g;\n',varname(FREQUENCY_LOW),FREQUENCY_LOW);
fprintf(fid,'%s=%g;\n',varname(FREQUENCY_HIGH),FREQUENCY_HIGH);
fprintf(fid,'%s=[%g %g];\n',varname(CONVOLUTION_SIZE),CONVOLUTION_SIZE);
fprintf(fid,'%s=%g;\n',varname(MINIMUM_OBJECT_AREA),MINIMUM_OBJECT_AREA);
fprintf(fid,'%s=%g;\n',varname(MERGE_HARMONICS),MERGE_HARMONICS);
fprintf(fid,'%s=%g;\n',varname(MERGE_HARMONICS_OVERLAP),MERGE_HARMONICS_OVERLAP);
fprintf(fid,'%s=%g;\n',varname(MERGE_HARMONICS_RATIO),MERGE_HARMONICS_RATIO);
fprintf(fid,'%s=%g;\n',varname(MERGE_HARMONICS_FRACTION),MERGE_HARMONICS_FRACTION);
fprintf(fid,'%s=%g;\n',varname(MINIMUM_VOCALIZATION_LENGTH),MINIMUM_VOCALIZATION_LENGTH);
fprintf(fid,'%s=[%s];\n',varname(CHANNELS),num2str(CHANNELS));
fclose(fid);



function ax2_print(i,left,right,type,filename,FS,NFFT,SAVE_WAV,SAVE_PNG,b,a,CHANNELS)

if isempty(CHANNELS)
  d=dir([filename '.ch*']);
  [tmp{1:length(d)}]=deal(d.name);
  CHANNELS=cellfun(@(x) str2num(x(end)),tmp);
end

tmp=[];  p=[];
for j=CHANNELS
  fid=fopen([filename '.ch' num2str(j)],'r');
  if(fid<0)  continue;  end
  fseek(fid,round(4*(round((left-0.025)*FS))),-1);
  tmp=fread(fid,round((right-left+0.050)*FS),'float32');
  fclose(fid);
  if(SAVE_WAV)
    tmp2=filtfilt(b,a,tmp);
    tmp2=tmp2./max([max(tmp2)-min(tmp2)]);
    wavwrite(tmp2,22000,[type num2str(i) '.ch' num2str(j) '.wav']);
  end
  [s,f,t,p(j,:,:)]=spectrogram(tmp,NFFT,[],[],FS,'yaxis');
end
tmp=squeeze(max(p,[],1));
tmp=log10(abs(tmp));
tmp4=reshape(tmp,1,prod(size(tmp)));
tmp2=prctile(tmp4,1);
tmp3=prctile(tmp4,99);
idx=find(tmp<tmp2);  tmp(idx)=tmp2;
idx=find(tmp>tmp3);  tmp(idx)=tmp3;
h=surf(t+left-0.025,f-f(2)/2,tmp,'EdgeColor','none');
uistack(h,'bottom');
colormap(gray);
axis tight;
set(gca,'xlim',[left right]+[-0.025 0.025]);
title([type ' #' num2str(i)]);
drawnow;
if(SAVE_PNG)  print('-dpng',[type num2str(i) '.png']);  end
delete(h);
