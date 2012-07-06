function mtbp2(filepath)

tmp=dir(fullfile(filepath,'*MTBP*.mat'));
idx=strfind({tmp.name},'_MTBP');
for i=1:length(tmp)
  datafiles{i}=tmp(i).name(1:(idx{i}(end)-1));
end
datafiles=unique(datafiles);
for i=1:length(datafiles)
  disp(fullfile(filepath,datafiles{i}));
  mtbp2_guts(fullfile(filepath,datafiles{i}));
end


function mtbp2_guts(filename)

if(exist('matlabpool')==2 && matlabpool('size')==0)
  matlabpool open
end

CHUNK_LEN=10;  % sec
OBJ_SIZE=1500;  % pixels
CONV_SIZE=[15 7];  % pixels
OVERLAP=0.005;  % sec

GROUNDTRUTH=0;
SAVE_PNG=0;
SAVE_WAV=0;

endings={'_MTBP128.mat','_MTBP256.mat','_MTBP512.mat'};


if(1)
disp('loading data...');
for i=1:length(endings)
  data(i)=load([filename endings{i}]);
end
if(GROUNDTRUTH)
groundtruth=load('../groundtruth.txt');
groundtruth=groundtruth(:,2:3);
groundtruth=sortrows(groundtruth,1);
idx=find(groundtruth(:,1)>groundtruth(:,2));
if(~isempty(idx))
  error(['rows ' num2str(idx) ' of groundtruth are invalid.']);
end
end
end


if(1)
disp('collapsing across channels and window sizes...');
df=data(3).df/10;
maxF=round(max([data(1).MT(:,2);  data(2).MT(:,2);  data(3).MT(:,2)]./df))+2+3+3;
maxT=max([data(1).MT(:,1);  2*data(2).MT(:,1)+1;  4*data(3).MT(:,1)+3]);

syls_chunked={};
chunk_len_win=round(CHUNK_LEN*data(1).FS/(data(1).NFFT/2));
parfor j=1:floor(maxT/chunk_len_win)
  chunk_idx=1+(j-1)*chunk_len_win;
%  disp(chunk_idx);

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

  im=im((round(20e3/df):round(120e3/df))+3+3,:);
  im=logical(conv2(single(im),ones(CONV_SIZE),'same'));

%% hough transform
%rotI=im(780:1050,188300:188490);
%%BW = edge(rotI,'canny');
%[H,theta,rho] = hough(rotI);
%P = houghpeaks(H,100,'threshold',ceil(0.2*max(H(:))),'nhoodsize',8*[15 5]+1);
%lines = houghlines(rotI,theta,rho,P,'FillGap',5,'MinLength',2*7);
%
%figure, imshow(rotI), hold on
%max_len = 0;
%for k = 1:length(lines)
%   xy = [lines(k).point1; lines(k).point2];
%   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%
%   % Plot beginnings and ends of lines
%   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%end

  syls_chunked{j}=bwconncomp(im,8);
end


disp('segmenting syllables...');
syls_separate.NumObjects=0;
syls_separate.PixelIdxList={};
syls_separate.Connectivity=syls_chunked{1}.Connectivity;
syls_separate.ImageSize=syls_chunked{1}.ImageSize;
for i=1:length(syls_chunked)
  syls_separate.NumObjects=syls_separate.NumObjects+syls_chunked{i}.NumObjects;
  for j=1:length(syls_chunked{i}.PixelIdxList)
    syls_chunked{i}.PixelIdxList{j}=syls_chunked{i}.PixelIdxList{j}+(i-1)*chunk_len_win*syls_chunked{i}.ImageSize(1);
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
skytruth_separate=[[tmp(:,1) tmp(:,1)+tmp(:,3)]*data(1).NFFT/2 zeros(size(tmp,1),1) [tmp(:,2) tmp(:,2)+tmp(:,4)].*df+20e3];


disp('calculating frequency contours...');
freq_contour={};
parfor i=1:length(syls_separate2)
  idx={};
  for j=1:3
    idx{j}=find(((2^(j-1)*data(j).MT(:,1))>=syls_separate2(i).BoundingBox(1)) & ...
                ((2^(j-1)*data(j).MT(:,1))<=sum(syls_separate2(i).BoundingBox([1 3]))) & ...
                (data(j).MT(:,2)>=(syls_separate2(i).BoundingBox(2)*df+20e3)) & ...
                (data(j).MT(:,2)<=(sum(syls_separate2(i).BoundingBox([2 4]))*df+20e3)));
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

end


if(1)
disp('merging nearby syllables...');
syls_merged=syls_separate;
syls_merged2=syls_separate2;
OVERLAP=OVERLAP*data(1).FS/data(1).NFFT*2;
tic;
for i=1:(length(syls_merged2)-1)
  if(isempty(syls_merged.PixelIdxList{i}))  continue;  end
  if(toc>10)  disp([num2str(i) ' of ' num2str(length(syls_merged2))]);  tic;  end;
  for j=(i+1):length(syls_merged2)
    if(isempty(syls_merged.PixelIdxList{j}))  continue;  end
    if(((sum(syls_merged2(i).BoundingBox([1 3]))+OVERLAP) > syls_merged2(j).BoundingBox(1)) &&...
       ((sum(syls_merged2(j).BoundingBox([1 3]))+OVERLAP) > syls_merged2(i).BoundingBox(1)))
      syls_merged.NumObjects=syls_merged.NumObjects-1;
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
    end
  end
end
idx=find(~cellfun(@isempty,syls_merged.PixelIdxList));
syls_merged.PixelIdxList={syls_merged.PixelIdxList{idx}};
syls_merged2=regionprops(syls_merged,'basic');
tmp=reshape([syls_merged2.BoundingBox],4,length(syls_merged2))';
skytruth_merged=[[tmp(:,1) tmp(:,1)+tmp(:,3)]*data(1).NFFT/2 zeros(size(tmp,1),1) [tmp(:,2) tmp(:,2)+tmp(:,4)].*df+20e3];
end


%skytruth=skytruth_separate;
skytruth=skytruth_merged;

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
disp([num2str(size(groundtruth,1)) ' manually segmented syllables, ' num2str(length(misses)) ' (' num2str(100*length(misses)/size(groundtruth,1),3) '%) of which are missed']);

false_alarms=find(skytruth(:,3)==0);
disp([num2str(size(skytruth,1)) ' automatically segmented syllables, ' num2str(length(false_alarms)) ' (' num2str(100*length(false_alarms)/size(skytruth,1),3) '%) of which are false alarms']);
end


%%%  plot
if(0)
disp('plotting...');
figure;
get(gcf,'position');
set(gcf,'position',[ans(1) ans(2) 4*ans(3) ans(4)]);

subplot('position',[0.05 0.1 0.9 0.8]);
set(gca,'color',[0 0 0]);
hold on;

%syls2=syls_separate2;
syls2=syls_merged2;
tmp=reshape([syls2.BoundingBox],4,length(syls2))';
tmp(:,1)=tmp(:,1)+(CONV_SIZE(2)-1)/2;
tmp(:,2)=tmp(:,2)+(CONV_SIZE(1)-1)/2;
tmp(:,3)=tmp(:,3)-CONV_SIZE(2);
tmp(:,4)=tmp(:,4)-CONV_SIZE(1);
plot([tmp(:,1) tmp(:,1)+tmp(:,3) tmp(:,1)+tmp(:,3) tmp(:,1) tmp(:,1)]'.*data(1).NFFT/2/data(1).FS,...
     [tmp(:,2) tmp(:,2) tmp(:,2)+tmp(:,4) tmp(:,2)+tmp(:,4) tmp(:,2)]'.*df+20e3,'y');

line([groundtruth(misses,1) groundtruth(misses,2) groundtruth(misses,2) groundtruth(misses,1) groundtruth(misses,1)]./data(1).FS,...
    [20e3 20e3 120e3 120e3 20e3],'color',[0 1 0]);
tmp=setdiff(1:size(groundtruth,1),misses);
line([groundtruth(tmp,1) groundtruth(tmp,2) groundtruth(tmp,2) groundtruth(tmp,1) groundtruth(tmp,1)]./data(1).FS,...
    [20e3 20e3 120e3 120e3 20e3],'color',[1 1 1]);

for i=1:length(freq_contour)
  plot(freq_contour{i}(:,1),freq_contour{i}(:,2),'r-');
end

axis tight;
xlabel('time (s)');
ylabel('frequency (Hz)');

if(0)  % plot hits, misses and false alarms separately
[b,a]=butter(4,20e3/(data(1).FS/2),'high');

hits=setdiff(1:size(groundtruth,1),misses);
for i=1:min([100 length(hits)])
  left =min([groundtruth(hits(i),1) skytruth(groundtruth(hits(i),3),1)]);
  right=max([groundtruth(hits(i),2) skytruth(groundtruth(hits(i),3),2)]);
  tmp=[];  p=[];
  for j=1:4
    fid=fopen(['../groundtruth.ch' num2str(j)],'r');
    fseek(fid,round(4*(left-round(0.025*data(1).FS))),-1);
    tmp=fread(fid,round(right-left+0.050*data(1).FS),'float32');
    fclose(fid);
    if(SAVE_WAV)
      tmp2=filtfilt(b,a,tmp);
      tmp2=tmp2./max([max(tmp2)-min(tmp2)]);
      wavwrite(tmp2,22000,['voclist' num2str(i) '.ch' num2str(j) '.wav']);
    end
    [s,f,t,p(j,:,:)]=spectrogram(tmp,data(2).NFFT,[],[],data(2).FS,'yaxis');
  end
  tmp=squeeze(max(p,[],1));
  tmp=log10(abs(tmp));
  tmp4=reshape(tmp,1,prod(size(tmp)));
  tmp2=prctile(tmp4,1);
  tmp3=prctile(tmp4,99);
  idx=find(tmp<tmp2);  tmp(idx)=tmp2;
  idx=find(tmp>tmp3);  tmp(idx)=tmp3;
  h=surf(t+left./data(1).FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
  colormap(gray);
  set(gca,'xlim',[left right]./data(1).FS+[-0.025 0.025]);
  drawnow;
  if(SAVE_PNG)  print('-dpng',['voclist' num2str(i) '.png']);  end
  delete(h);
end

for i=1:min([100 length(misses)])
  tmp=[];  p=[];
  for j=1:4
    fid=fopen(['../groundtruth.ch' num2str(j)],'r');
    fseek(fid,round(4*(groundtruth(misses(i),1)-round(0.025*data(1).FS))),-1);
    tmp=fread(fid,round(diff(groundtruth(misses(i),1:2))+0.050*data(1).FS),'float32');
    fclose(fid);
    if(SAVE_WAV)
      tmp2=filtfilt(b,a,tmp);
      tmp2=tmp2./max([max(tmp2) -min(tmp2)]);
      wavwrite(tmp2,22000,['miss' num2str(i) '.ch' num2str(j) '.wav']);
    end
    [s,f,t,p(j,:,:)]=spectrogram(tmp,data(2).NFFT,[],[],data(2).FS,'yaxis');
  end
  tmp=squeeze(max(p,[],1));
  tmp=log10(abs(tmp));
  tmp4=reshape(tmp,1,prod(size(tmp)));
  tmp2=prctile(tmp4,1);
  tmp3=prctile(tmp4,99);
  idx=find(tmp<tmp2);  tmp(idx)=tmp2;
  idx=find(tmp>tmp3);  tmp(idx)=tmp3;
  h=surf(t+groundtruth(misses(i),1)./data(1).FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
  colormap(gray);
  set(gca,'xlim',groundtruth(misses(i),1:2)./data(1).FS+[-0.025 0.025]);
  drawnow;
  if(SAVE_PNG)  print('-dpng',['miss' num2str(i) '.png']);  end
  delete(h);
end

for i=1:min([100 length(false_alarms)])
  tmp=[];  p=[];
  for j=1:4
    fid=fopen(['../groundtruth.ch' num2str(j)],'r');
    fseek(fid,round(4*(skytruth(false_alarms(i),1)-round(0.025*data(1).FS))),-1);
    tmp=fread(fid,round(diff(skytruth(false_alarms(i),1:2))+0.050*data(1).FS),'float32');
    fclose(fid);
    if(SAVE_WAV)
      tmp2=filtfilt(b,a,tmp);
      tmp2=tmp2./max([max(tmp2) -min(tmp2)]);
      wavwrite(tmp2,22000,['false_alarm' num2str(i) '.ch' num2str(j) '.wav']);
    end
    [s,f,t,p(j,:,:)]=spectrogram(tmp,data(2).NFFT,[],[],data(2).FS,'yaxis');
  end
  tmp=squeeze(max(p,[],1));
  tmp=log10(abs(tmp));
  tmp4=reshape(tmp,1,prod(size(tmp)));
  tmp2=prctile(tmp4,1);
  tmp3=prctile(tmp4,99);
  idx=find(tmp<tmp2);  tmp(idx)=tmp2;
  idx=find(tmp>tmp3);  tmp(idx)=tmp3;
  h=surf(t+skytruth(false_alarms(i),1)./data(1).FS-0.025,f-f(2)/2,tmp,'EdgeColor','none');
  colormap(gray);
  set(gca,'xlim',skytruth(false_alarms(i),1:2)./data(1).FS+[-0.025 0.025]);
  drawnow;
  if(SAVE_PNG)  print('-dpng',['false_alarm' num2str(i) '.png']);  end
  delete(h);
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


%fid=fopen('frqlist.txt','w');
%for i=1:size(skytruth,1)
%  j=skytruth(i,1);
%  while(j<=skytruth(i,2))
%    %fwrite(fid,'%f',data(1).MT{});
%  end
%  fwrite(fid,',');
%end
%fclose(fid);
