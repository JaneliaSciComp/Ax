#!/bin/bash

# ./convert.sh full_path_to_directory_of_-out_folders

dir_name=$1
data_names=$(ls -1 $dir_name | grep \\-out)
for d in $data_names
do
  data_prefix=${d%-out*}
  file_names=$(ls -1 $dir_name/$d)
  for f in $file_names
  do
    file_prefix=${f%.*}
    echo 'cp' $dir_name/$d/$f $dir_name/${data_prefix}.${file_prefix}1234
    cp $dir_name/$d/$f $dir_name/${data_prefix}.${file_prefix}1234
  done
done
