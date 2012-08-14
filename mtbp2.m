%function mtbp2(type,merge,filepath)
%
%type='R' is for rejection calls, 'US' for ultrasonic
%set merge to the overlap criterion, in seconds, below which vocalizations
%  are combined, or to [] to not combine
%
%mtbp2('US',[],'/groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/');
%mtbp2('R',0.005,'/groups/egnor/egnorlab/for_ben/sys_test_07052012a/demux/');

function mtbp2(type,merge,filepath)

if((nargin~=3) || ~ismember(upper(type),{'R' 'US'}))
  error('invalid args');
end

tmp=dir(fullfile(filepath,'*_MTBP*.mat'));
idx=strfind({tmp.name},'_MTBP');
for i=1:length(tmp)
  datafiles{i}=tmp(i).name(1:(idx{i}(end)-1));
end
datafiles=unique(datafiles);
for i=1:length(datafiles)
  disp(fullfile(filepath,datafiles{i}));
  mtbp2_guts(type,merge,fullfile(filepath,datafiles{i}));
end


function mtbp2_guts(type,merge,filename)

if(exist('matlabpool')==2 && matlabpool('size')==0)
  matlabpool open
end

CHUNK_LEN=10;  % sec
GROUNDTRUTH=0;
SAVE_PNG=1;
SAVE_WAV=0;
MERGE=merge;

switch(upper(type))
  case 'US'   % ultrasonic
    OBJ_SIZE=1500;  % pixels
    CONV_SIZE=[15 7];  % pixels, must be odd, should convert this to Hz x sec
    F_LOW=20e3;  % Hz
    F_HIGH=120e3;  % Hz
    NHARM=1;

  case 'R'   % rejection
    OBJ_SIZE=1000;  % pixels
    CONV_SIZE=[15 7];  % pixels, must be odd, should convert this to Hz x sec
    F_LOW=1e3;  % Hz
    F_HIGH=20e3;  % Hz
    NHARM=3;
end


disp('loading data...');
tmp=dir([filename,'_MTBP*.mat']);
for i=1:length(tmp)
  data(i)=load(fullfile(fileparts(filename),tmp(i).name));
end
if(GROUNDTRUTH)
  groundtruth=load([filename '.txt']);
  groundtruth=groundtruth(:,2:3);
  groundtruth=sortrows(groundtruth,1);
  idx=find(groundtruth(:,1)>groundtruth(:,2));
  if(~isempty(idx))
    error(['rows ' num2str(idx) ' of groundtruth are invalid.']);
  end
end


disp('collapsing across channels and window sizes...');
df=data(3).df/10;
maxF=round(max([data(1).MT(:,2);  data(2).MT(:,2);  data(3).MT(:,2)]./df))+2+3+3;
maxT=max([data(1).MT(:,1);  2*data(2).MT(:,1)+1;  4*data(3).MT(:,1)+3]);

syls_chunked={};
chunk_len=round(CHUNK_LEN*data(1).FS/(data(1).NFFT/2));

maxT/round(maxT/chunk_len);
chunk_splits=round(1:ans:(maxT+1));
parfor j=1:(length(chunk_splits)-1)
  chunk_idx=chunk_splits(j);
  chunk_len_win=chunk_splits(j+1)-chunk_splits(j)-1;

  im=false(maxF,chunk_len_win);
  idx=find((data(1).MT(:,1)>=chunk_idx) & (data(1).MT(:,1)<(chunk_idx+chunk_len_win)));
  for i=((-2:2)+3+3)
    im(sub2ind([maxF,chunk_len_win],round(data(1).MT(idx,2)/df  )+i,  data(1).MT(idx,1)-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(1).MT(idx,2)/df-1)+i,  data(1).MT(idx,1)-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(1).MT(idx,2)/df-2)+i,  data(1).MT(idx,1)-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(1).MT(idx,2)/df-3)+i,  data(1).MT(idx,1)-chunk_idx+1))=true;
  end
  idx=find(((2*data(2).MT(:,1)-1)>=chunk_idx) & ((2*data(2).MT(:,1)+1)<(chunk_idx+chunk_len_win)));
  for i=((-1:1)+3+3)
    im(sub2ind([maxF,chunk_len_win],round(data(2).MT(idx,2)/df  )+i,2*data(2).MT(idx,1)-1-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(2).MT(idx,2)/df  )+i,2*data(2).MT(idx,1)  -chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(2).MT(idx,2)/df  )+i,2*data(2).MT(idx,1)+1-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(2).MT(idx,2)/df-1)+i,2*data(2).MT(idx,1)-1-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(2).MT(idx,2)/df-1)+i,2*data(2).MT(idx,1)  -chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(2).MT(idx,2)/df-1)+i,2*data(2).MT(idx,1)+1-chunk_idx+1))=true;
  end
  idx=find(((4*data(3).MT(:,1)-3)>=chunk_idx) & ((4*data(3).MT(:,1)+3)<(chunk_idx+chunk_len_win)));
  for i=0
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)-3-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)-2-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)-1-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)  -chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)+1-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)+2-chunk_idx+1))=true;
    im(sub2ind([maxF,chunk_len_win],round(data(3).MT(idx,2)/df  )+i,4*data(3).MT(idx,1)+3-chunk_idx+1))=true;
  end

  tmp=(round(F_LOW/df):round(F_HIGH/df))+3+3;
  im=[zeros(length(tmp),(CONV_SIZE(2)-1)/2) im(tmp,:) zeros(length(tmp),(CONV_SIZE(2)-1)/2)];
  im=logical(conv2(single(im),ones(CONV_SIZE),'same'));

  syls_chunked{j}=bwconncomp(im,8);
end


disp('segmenting syllables...');
tic;
for k=1:(length(syls_chunked)-1)
  flag=1;
  while(flag)
    flag=0;
    for i=1:syls_chunked{k}.NumObjects
      if(toc>10)  disp([num2str(k) ' of ' num2str(length(syls_chunked))]);  tic;  end;
      [ri ci]=ind2sub(syls_chunked{k}.ImageSize,syls_chunked{k}.PixelIdxList{i});
      if(sum(ci>=(syls_chunked{k}.ImageSize(2)-(CONV_SIZE(2)-1)/2))==0) continue;  end
      j=1;
      while j<=syls_chunked{k+1}.NumObjects
        [rj cj]=ind2sub(syls_chunked{k+1}.ImageSize,syls_chunked{k+1}.PixelIdxList{j});
        if(sum(cj<=((CONV_SIZE(2)-1)/2))==0)  j=j+1;  continue;  end
        %if((max(ri) < (min(rj)-1)) || (max(rj) < (min(ri)-1)))  j=j+1;  continue;  end
        cj=cj+chunk_splits(k+1)-chunk_splits(k);
        min(min((repmat(ri,1,length(rj))-repmat(rj',length(ri),1)).^2 + ...
                (repmat(ci,1,length(cj))-repmat(cj',length(ci),1)).^2));
        if ans<=2
          disp(['unsplitting syllable between chunks #' num2str(k) '-' num2str(k+1)]);
          flag=1;
          syls_chunked{k}.PixelIdxList{i}=[syls_chunked{k}.PixelIdxList{i}; ...
              syls_chunked{k+1}.PixelIdxList{j}+(chunk_splits(k+1)-chunk_splits(k))*syls_chunked{k}.ImageSize(1)];
          syls_chunked{k+1}.PixelIdxList(j)=[];
          syls_chunked{k+1}.NumObjects=syls_chunked{k+1}.NumObjects-1;
        else
          j=j+1;
        end
      end
    end
  end
end

syls_separate.NumObjects=0;
syls_separate.PixelIdxList={};
syls_separate.Connectivity=syls_chunked{1}.Connectivity;
syls_separate.ImageSize=syls_chunked{1}.ImageSize;
for i=1:length(syls_chunked)
  syls_separate.NumObjects=syls_separate.NumObjects+syls_chunked{i}.NumObjects;
  for j=1:length(syls_chunked{i}.PixelIdxList)
    %syls_chunked{i}.PixelIdxList{j}=syls_chunked{i}.PixelIdxList{j}+(i-1)*chunk_len_win*syls_chunked{i}.ImageSize(1);
    syls_chunked{i}.PixelIdxList{j}=syls_chunked{i}.PixelIdxList{j}+(chunk_splits(i)-1)*syls_chunked{i}.ImageSize(1);
  end
  syls_separate.PixelIdxList=[syls_separate.PixelIdxList syls_chunked{i}.PixelIdxList];
end

syls_separate2=regionprops(syls_separate,'basic');
tmp=find([syls_separate2.Area]>OBJ_SIZE);
syls_separate2=syls_separate2(tmp);
syls_separate.PixelIdxList=syls_separate.PixelIdxList(tmp);
syls_separate.NumObjects=length(tmp);
tmp=reshape([syls_separate2.BoundingBox],4,length(syls_separate2))';
tmp(:,1)=tmp(:,1)+(CONV_SIZE(2)-1)/2;
tmp(:,2)=tmp(:,2)+(CONV_SIZE(1)-1)/2;
tmp(:,3)=tmp(:,3)-CONV_SIZE(2);
tmp(:,4)=tmp(:,4)-CONV_SIZE(1);
skytruth_separate=[[tmp(:,1) tmp(:,1)+tmp(:,3)]*data(1).NFFT/2 zeros(size(tmp,1),1) [tmp(:,2) tmp(:,2)+tmp(:,4)].*df+F_LOW];


disp('calculating frequency contours...');
freq_contour={};
parfor i=1:length(syls_separate2)
  idx={};
  for j=1:3
    idx{j}=find(((2^(j-1)*data(j).MT(:,1))>=syls_separate2(i).BoundingBox(1)) & ...
                ((2^(j-1)*data(j).MT(:,1))<=sum(syls_separate2(i).BoundingBox([1 3]))) & ...
                (data(j).MT(:,2)>=(syls_separate2(i).BoundingBox(2)*df+F_LOW)) & ...
                (data(j).MT(:,2)<=(sum(syls_separate2(i).BoundingBox([2 4]))*df+F_LOW)));
  end
  tmp=[bsxfun(@times,data(1).MT(idx{1},1:3),[data(1).NFFT/2/data(1).FS 1 1]); ...
       bsxfun(@times,data(2).MT(idx{2},1:3),[data(2).NFFT/2/data(1).FS 1 1]);...
       bsxfun(@times,data(3).MT(idx{3},1:3),[data(3).NFFT/2/data(1).FS 1 1])];
  tmp=sortrows(tmp);
  freq_contour{i}=zeros(length(unique(tmp(:,1))),3);
  j=1;  l=1;
  while(j<=size(tmp,1))
    k=j+1;  while((k<=size(tmp,1)) && (tmp(j,1)==tmp(k,1)))  k=k+1;  end
    [~, idx]=max(tmp(j:(k-1),3));
    freq_contour{i}(l,:)=tmp(j+idx-1,:);
    j=k;  l=l+1;
  end
end


if(~isempty(MERGE))
  disp('merging nearby syllables...');
  syls_merged=syls_separate;
  syls_merged2=syls_separate2;
  syls_merged3=ones(1,length(syls_merged2));
  MERGE=MERGE*data(1).FS/data(1).NFFT*2;
  tic;
  for i=1:(length(syls_merged2)-1)
    if(isempty(syls_merged.PixelIdxList{i}))  continue;  end
    if(toc>10)  disp([num2str(i) ' of ' num2str(length(syls_merged2))]);  tic;  end;
    flag=1;
    while(flag)
      flag=0;
      for j=(i+1):length(syls_merged2)
        if(isempty(syls_merged.PixelIdxList{j}))  continue;  end
        if(((sum(syls_merged2(i).BoundingBox([1 3]))+MERGE) > syls_merged2(j).BoundingBox(1)) &&...
           ((sum(syls_merged2(j).BoundingBox([1 3]))+MERGE) > syls_merged2(i).BoundingBox(1)))
          flag=1;
%          syls_merged.NumObjects=syls_merged.NumObjects-1;
          syls_merged.PixelIdxList{i}=[syls_merged.PixelIdxList{i}; syls_merged.PixelIdxList{j}];
          syls_merged.PixelIdxList{j}=[];
          syls_merged2(i).BoundingBox(1)=min([syls_merged2(i).BoundingBox(1) syls_merged2(j).BoundingBox(1)]);
          syls_merged2(i).BoundingBox(2)=min([syls_merged2(i).BoundingBox(2) syls_merged2(j).BoundingBox(2)]);
          syls_merged2(i).BoundingBox(3)=...
              max([sum(syls_merged2(i).BoundingBox([1 3])) sum(syls_merged2(j).BoundingBox([1 3]))])-...
              syls_merged2(i).BoundingBox(1);
          syls_merged2(i).BoundingBox(4)=...
              max([sum(syls_merged2(i).BoundingBox([2 4])) sum(syls_merged2(j).BoundingBox([2 4]))])-...
              syls_merged2(i).BoundingBox(2);
          syls_merged3(i)=syls_merged3(i)+syls_merged3(j);
          syls_merged3(j)=0;
        end
      end
    end
  end
  %idx=find(~cellfun(@isempty,syls_merged.PixelIdxList));
  idx=find(syls_merged3>=NHARM);
  syls_merged.NumObjects=length(idx);
  syls_merged.PixelIdxList={syls_merged.PixelIdxList{idx}};
  syls_merged2=regionprops(syls_merged,'basic');
  tmp=reshape([syls_merged2.BoundingBox],4,length(syls_merged2))';
  skytruth=[[tmp(:,1) tmp(:,1)+tmp(:,3)]*data(1).NFFT/2 zeros(size(tmp,1),1) [tmp(:,2) tmp(:,2)+tmp(:,4)].*df+F_LOW];
  syls2=syls_merged2;
else
  tmp=reshape([syls_separate2.BoundingBox],4,length(syls_separate2))';
  skytruth=[[tmp(:,1) tmp(:,1)+tmp(:,3)]*data(1).NFFT/2 zeros(size(tmp,1),1) [tmp(:,2) tmp(:,2)+tmp(:,4)].*df+F_LOW];
  syls2=syls_separate2;
end


if(GROUNDTRUTH)
  disp('comparing to ground truth...');  drawnow;
  groundtruth(:,3)=0;
  skytruth(:,3)=0;
  a=1;  m=1;
  while((a<size(skytruth,1))&&(skytruth(a,2)<groundtruth(m,1)))
    a=a+1;
  end
  %tic;
  while((a<=size(skytruth,1)) && (m<=size(groundtruth,1)))
  %  if(toc>1)  disp([num2str(a) ' ' num2str(m)]);  tic;  end
    if(skytruth(a,1)<groundtruth(m,2))
      groundtruth(m,3)=a;
      skytruth(a,3)=m;
      a=a+1;
    end
    if(a>size(skytruth,1))  break;  end
    if(skytruth(a,1)>groundtruth(m,2))
      m=m+1;
      if(m>size(groundtruth,1))  break;  end
    end
    while((a<size(skytruth,1))&&(skytruth(a,2)<groundtruth(m,1)))
      a=a+1;
    end
  end

  misses=find(groundtruth(:,3)==0);
  disp([num2str(size(groundtruth,1)) ' manually segmented syllables, ' num2str(length(misses)) ...
      ' (' num2str(100*length(misses)/size(groundtruth,1),3) '%) of which are missed']);

  false_alarms=find(skytruth(:,3)==0);
  disp([num2str(size(skytruth,1)) ' automatically segmented syllables, ' num2str(length(false_alarms)) ...
      ' (' num2str(100*length(false_alarms)/size(skytruth,1),3) '%) of which are false alarms']);
else
  disp([num2str(size(skytruth,1)) ' automatically segmented syllables']);
end


if(SAVE_WAV || SAVE_PNG)
disp('plotting...');
figure;
get(gcf,'position');
set(gcf,'position',[ans(1) ans(2) 4*ans(3) ans(4)]);

subplot('position',[0.05 0.1 0.9 0.8]);
set(gca,'color',[0 0 0]);
hold on;

tmp=reshape([syls2.BoundingBox],4,length(syls2))';
tmp(:,1)=tmp(:,1)+(CONV_SIZE(2)-1)/2;
tmp(:,2)=tmp(:,2)+(CONV_SIZE(1)-1)/2;
tmp(:,3)=tmp(:,3)-CONV_SIZE(2);
tmp(:,4)=tmp(:,4)-CONV_SIZE(1);
plot([tmp(:,1) tmp(:,1)+tmp(:,3) tmp(:,1)+tmp(:,3) tmp(:,1) tmp(:,1)]'.*data(1).NFFT/2/data(1).FS,...
     [tmp(:,2) tmp(:,2) tmp(:,2)+tmp(:,4) tmp(:,2)+tmp(:,4) tmp(:,2)]'.*df+F_LOW,'y');

if(GROUNDTRUTH)
line([groundtruth(misses,1) groundtruth(misses,2) groundtruth(misses,2) groundtruth(misses,1) groundtruth(misses,1)]./data(1).FS,...
    [F_LOW F_LOW F_HIGH F_HIGH F_LOW],'color',[0 1 0]);
tmp=setdiff(1:size(groundtruth,1),misses);
line([groundtruth(tmp,1) groundtruth(tmp,2) groundtruth(tmp,2) groundtruth(tmp,1) groundtruth(tmp,1)]./data(1).FS,...
    [F_LOW F_LOW F_HIGH F_HIGH F_LOW],'color',[1 1 1]);
end

for i=1:length(freq_contour)
  plot(freq_contour{i}(:,1),freq_contour{i}(:,2),'r-');
end

plot(repmat(chunk_splits*data(1).NFFT/2/data(1).FS,2,1),...
     repmat([F_LOW; F_HIGH],1,length(chunk_splits)),'c');

axis tight;
xlabel('time (s)');
ylabel('frequency (Hz)');

[b,a]=butter(4,F_LOW/(data(1).FS/2),'high');


% plot hits, misses and false alarms separately
if(~GROUNDTRUTH)
  for i=1:min([20 length(skytruth)])
    left=skytruth(i,1);
    right=skytruth(i,2);
    mtbp2_print(i,left,right,'voc',filename,data(1).FS,data(1).NFFT,SAVE_WAV,SAVE_PNG);
  end
else

hits=setdiff(1:size(groundtruth,1),misses);
for i=1:min([20 length(hits)])
  left =min([groundtruth(hits(i),1) skytruth(groundtruth(hits(i),3),1)]);
  right=max([groundtruth(hits(i),2) skytruth(groundtruth(hits(i),3),2)]);
  mtbp2_print(i,left,right,'hit',filename,data(1).FS,data(1).NFFT,SAVE_WAV,SAVE_PNG);

%  tmp=[];  p=[];
%  for j=1:4
%    fid=fopen([filename '.ch' num2str(j)],'r');
%    fseek(fid,round(4*(left-round(0.025*data(1).FS))),-1);
%    tmp=fread(fid,round(right-left+0.050*data(1).FS),'float32');
%    fclose(fid);
%    if(SAVE_WAV)
%      tmp2=filtfilt(b,a,tmp);
%      tmp2=tmp2./max([max(tmp2)-min(tmp2)]);
%      wavwrite(tmp2,22000,[fileparts(filename) '/voclist' num2str(i) '.ch' num2str(j) '.wav']);
%    end
%    [s,f,t,p(j,:,:)]=spectrogram(tmp,data(2).NFFT,[],[],data(2).FS,'yaxis');
%  end
%  tmp=squeeze(max(p,[],1));
%  tmp=log10(abs(tmp));
%  tmp4=reshape(tmp,1,prod(size(tmp)));
%  tmp2=prctile(tmp4,1);
%  tmp3=prctile(tmp4,99);
%  idx=find(tmp<tmp2);  tmp(idx)=tmp2;
%  idx=find(tmp>tmp3);  tmp(idx)=tmp3;
%  h=surf(t+left./data(1).FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
%  colormap(gray);
%  set(gca,'xlim',[left right]./data(1).FS+[-0.025 0.025]);
%  title(['hit #' num2str(i)]);
%  drawnow;
%  if(SAVE_PNG)  print('-dpng',[fileparts(filename) '/voclist' num2str(i) '.png']);  end
%  delete(h);
end

for i=1:min([100 length(misses)])
  left=groundtruth(misses(i),1);
  right=groundtruth(misses(i),2);
  mtbp2_print(i,left,right,'miss',filename,data(1).FS,data(1).NFFT,SAVE_WAV,SAVE_PNG);

%  tmp=[];  p=[];
%  for j=1:4
%    fid=fopen([filename '.ch' num2str(j)],'r');
%    fseek(fid,round(4*(groundtruth(misses(i),1)-round(0.025*data(1).FS))),-1);
%    tmp=fread(fid,round(diff(groundtruth(misses(i),1:2))+0.050*data(1).FS),'float32');
%    fclose(fid);
%    if(SAVE_WAV)
%      tmp2=filtfilt(b,a,tmp);
%      tmp2=tmp2./max([max(tmp2) -min(tmp2)]);
%      wavwrite(tmp2,22000,[fileparts(filename) '/miss' num2str(i) '.ch' num2str(j) '.wav']);
%    end
%    [s,f,t,p(j,:,:)]=spectrogram(tmp,data(2).NFFT,[],[],data(2).FS,'yaxis');
%  end
%  tmp=squeeze(max(p,[],1));
%  tmp=log10(abs(tmp));
%  tmp4=reshape(tmp,1,prod(size(tmp)));
%  tmp2=prctile(tmp4,1);
%  tmp3=prctile(tmp4,99);
%  idx=find(tmp<tmp2);  tmp(idx)=tmp2;
%  idx=find(tmp>tmp3);  tmp(idx)=tmp3;
%  h=surf(t+groundtruth(misses(i),1)./data(1).FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
%  colormap(gray);
%  set(gca,'xlim',groundtruth(misses(i),1:2)./data(1).FS+[-0.025 0.025]);
%  title(['miss #' num2str(i)]);
%  drawnow;
%  if(SAVE_PNG)  print('-dpng',[fileparts(filename) '/miss' num2str(i) '.png']);  end
%  delete(h);
end

for i=1:min([100 length(false_alarms)])
  left=skytruth(false_alarms(i),1);
  right=skytruth(false_alarms(i),2);
  mtbp2_print(i,left,right,'false_alarm',filename,data(1).FS,data(1).NFFT,SAVE_WAV,SAVE_PNG);

%  tmp=[];  p=[];
%  for j=1:4
%    fid=fopen([filename '.ch' num2str(j)],'r');
%    fseek(fid,round(4*(skytruth(false_alarms(i),1)-round(0.025*data(1).FS))),-1);
%    tmp=fread(fid,round(diff(skytruth(false_alarms(i),1:2))+0.050*data(1).FS),'float32');
%    fclose(fid);
%    if(SAVE_WAV)
%      tmp2=filtfilt(b,a,tmp);
%      tmp2=tmp2./max([max(tmp2) -min(tmp2)]);
%      wavwrite(tmp2,22000,[fileparts(filename) '/false_alarm' num2str(i) '.ch' num2str(j) '.wav']);
%    end
%    [s,f,t,p(j,:,:)]=spectrogram(tmp,data(2).NFFT,[],[],data(2).FS,'yaxis');
%  end
%  tmp=squeeze(max(p,[],1));
%  tmp=log10(abs(tmp));
%  tmp4=reshape(tmp,1,prod(size(tmp)));
%  tmp2=prctile(tmp4,1);
%  tmp3=prctile(tmp4,99);
%  idx=find(tmp<tmp2);  tmp(idx)=tmp2;
%  idx=find(tmp>tmp3);  tmp(idx)=tmp3;
%  h=surf(t+skytruth(false_alarms(i),1)./data(1).FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
%  colormap(gray);
%  set(gca,'xlim',skytruth(false_alarms(i),1:2)./data(1).FS+[-0.025 0.025]);
%  title(['false_alarm #' num2str(i)]);
%  drawnow;
%  if(SAVE_PNG)  print('-dpng',[fileparts(filename) '/false_alarm' num2str(i) '.png']);  end
%  delete(h);
end
end
end


disp('dumping files...');
tmp=[skytruth(:,1:2)./data(1).FS skytruth(:,4:5)];
save([filename '_voclist.txt'],'tmp','-ascii');
save([filename '_fcontours.mat'],'freq_contour');
if(GROUNDTRUTH)
tmp=groundtruth(misses,1:2)./data(1).FS;
save([filename '_misses.txt'],'tmp','-ascii');
tmp=[skytruth(false_alarms,1:2)./data(1).FS skytruth(false_alarms,4:5)];
save([filename '_false_alarms.txt'],'tmp','-ascii');
end



function mtbp2_print(i,left,right,type,filename,FS,NFFT,SAVE_WAV,SAVE_PNG)

tmp=[];  p=[];
for j=1:4
  fid=fopen([filename '.ch' num2str(j)],'r');
  fseek(fid,round(4*(left-round(0.025*FS))),-1);
  tmp=fread(fid,round(right-left+0.050*FS),'float32');
  fclose(fid);
  if(SAVE_WAV)
    tmp2=filtfilt(b,a,tmp);
    tmp2=tmp2./max([max(tmp2)-min(tmp2)]);
    wavwrite(tmp2,22000,[fileparts(filename) '/' type num2str(i) '.ch' num2str(j) '.wav']);
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
h=surf(t+left./FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
colormap(gray);
set(gca,'xlim',[left right]./FS+[-0.025 0.025]);
title([type ' #' num2str(i)]);
drawnow;
if(SAVE_PNG)  print('-dpng',[fileparts(filename) '/' type num2str(i) '.png']);  end
delete(h);
