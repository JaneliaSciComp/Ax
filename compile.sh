#ssh login
#qlogin -l interactive=true,matlab=1
#cd to ax/
#./compile.sh

#hard-coded for matlab 2012b on janelia cluster

mkdir -p ax

/usr/local/matlab-2012b/bin/mcc -o ax \
  -W main:ax \
  -T link:exe \
  -d ax \
  -w enable:specified_file_mismatch \
  -w enable:repeated_file \
  -w enable:switch_ignored \
  -w enable:missing_lib_sentinel \
  -w enable:demo_license \
  -v ax.m \
  -a chronux
