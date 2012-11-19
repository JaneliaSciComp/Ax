FS=450450;
NFFT=[0.009 0.0045 0.0022];
NW=18;
K=24;
PVAL=0.01;

channels=1:4;
obj_size=1000;
conv_size=[15 7];
f_low=1e3;
f_high=20e3;
nseg=3;
merge_time=[];
merge_freq=0;
  merge_freq_overlap=0.9;
  merge_freq_ratio=0.1;
  merge_freq_fraction=0.9;
