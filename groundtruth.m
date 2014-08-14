% function [miss, false_alarm]=groundtruth(ground_file, sky_file, time_range)
%
% ground_file and sky_file are voc list text files.  each row is a
% vocalization, the four tab-delimited columns are the start and stop times
% and low and high frequencies
%
% time_range is 2-element array specifing the subset of the files to
% consider in seconds.  use [] for the whole file.
%
% miss.area, miss.freq, and miss.time are vectors indicating the fraction
% of the area, frequency, and time of each manually annotated
% vocalizations' bounding box that is overlapped by automatically annotated
% vocalizations, and miss.num is the number of the latter which overlap the
% former (splits).
%
% false_alarm.area, false_alarm.freq, and false_alarm.time are vectors
% indicating the fraction of the area, frequency, and time of each
% automatically annotated vocalization's bounding box that is overlapped by
% manually annotated vocalizations, and false_alarm.num is the number of
% the latter which overlap the former (lumps).
%
% e.g.
% [miss, false_alarm] = groundtruth(...
% '/groups/egnor/egnorlab/ben/groundtruth/kelly_vs_meghan/m53606_m50733_urine_together_20121218_KellyAxGT20140311.voc',...
% '/groups/egnor/egnorlab/ben/groundtruth/kelly_vs_meghan/m53606_m50733_urine_together_20121218_MeganAxGT20140311.voc',...
% [90 105]);

function [miss, false_alarm]=groundtruth(ground_file, sky_file, time_range)

plot_flag=true;

ground_box=load(ground_file);  % a .txt file from e.g. Tempo
sky_box=load(sky_file);  % a .txt file from e.g. Ax

if(~isempty(time_range))
  idx=find((ground_box(:,1)>time_range(1)) & (ground_box(:,1)<time_range(2)) & (ground_box(:,2)>time_range(1)) & (ground_box(:,2)<time_range(2)));
  ground_box=ground_box(idx,:);
  if ~isempty(sky_box)
    idx=find((sky_box(:,1)>time_range(1)) & (sky_box(:,1)<time_range(2)) & (sky_box(:,2)>time_range(1)) & (sky_box(:,2)<time_range(2)));
    sky_box=sky_box(idx,:);
  end
end

if(plot_flag)
  figure;
  subplot(2,1,1);  hold on;
end
miss=groundtruth_guts(ground_box, sky_box, 1, plot_flag);
false_alarm=groundtruth_guts(sky_box, ground_box, 3, plot_flag);

if(plot_flag)
  subplot(2,2,3);
  [n,x]=hist(miss.area,50);
  bar(x,n,'barwidth',1,'facecolor','k');
  num=sum(miss.area>0.5);  den=length(miss.area);
  title([num2str(num) '/' num2str(den) '=' num2str(num/den*100,3) '% misses, ' ...
      num2str(sum(miss.num>1)/den*100,3) '% splits']);
  xlabel('miss score')
  ylabel('# vocs');
  axis tight;

  subplot(2,2,4);
  [n,x]=hist(false_alarm.area,50);
  bar(x,n,'barwidth',1,'facecolor','k');
  num=sum(false_alarm.area>0.5);  den=length(false_alarm.area);
  title([num2str(num) '/' num2str(den) '=' num2str(num/den*100,3) '% false alarms, ' ...
      num2str(sum(false_alarm.num>1)/den*100,3) '% lumps']);
  xlabel('false alarm score');
  ylabel('# vocs');
  axis tight;
end


function [x,y]=bbox2polygon(bbox)

x=[bbox(1) bbox(1) bbox(2) bbox(2) bbox(1)];
y=[bbox(3) bbox(4) bbox(4) bbox(3) bbox(3)];


function overlap = groundtruth_guts(this, that, color, plot_flag)

overlap.area=nan(size(this,1),1);
overlap.freq=nan(size(this,1),1);
overlap.time=nan(size(this,1),1);
overlap.num=zeros(size(this,1),1);
overlap.size=(this(:,2)-this(:,1)).*(this(:,4)-this(:,3));  % in Hz-sec
for idx_this=1:size(this,1)
  [x0,y0]=bbox2polygon(this(idx_this,:));
  x0p=x0;  y0p=y0;
  freqB=[]; freqT=[];
  timeL=[]; timeR=[];
  for idx_that=1:size(that,1)
    if ((this(idx_this,2)<that(idx_that,1)) || (this(idx_this,1)>that(idx_that,2)) || ...
        (this(idx_this,4)<that(idx_that,3)) || (this(idx_this,3)>that(idx_that,4)))  continue;  end
    [x,y]=bbox2polygon(that(idx_that,:));
    [x0p,y0p]=polybool('subtraction',x0p,y0p,x,y);
    [freqB freqT]=MergeBrackets([freqB that(idx_that,3)],[freqT that(idx_that,4)]);
    [timeL timeR]=MergeBrackets([timeL that(idx_that,1)],[timeR that(idx_that,2)]);
    overlap.num(idx_this) = overlap.num(idx_this) + 1;
  end
  [x0p,y0p]=polysplit(x0p,y0p);
  overlap.area(idx_this) = sum(cellfun(@(x,y) (ispolycw(x,y)*2-1)*polyarea(x,y), x0p,y0p)) / polyarea(x0,y0);
  [l r]=RangeIntersection(freqB,freqT,min(y0),max(y0));
  overlap.freq(idx_this) = ((max(y0)-min(y0)) - sum(r-l)) / (max(y0)-min(y0));
  [l r]=RangeIntersection(timeL,timeR,min(x0),max(x0));
  overlap.time(idx_this) = ((max(x0)-min(x0)) - sum(r-l)) / (max(x0)-min(x0));
  if(plot_flag)
    c=[0 0 0];  c(color)=max(0,overlap.area(idx_this));
    line(this(idx_this,[1 2 2 1 1]),this(idx_this,[3 3 4 4 3]),'color',c,'linewidth',2);
  end
end
