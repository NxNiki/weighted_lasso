#!/bin/sh
#set -e
current_dir=$(pwd)
time=$(date +%Y_%m_%d_%H_%M_%S)

backup_dir="/home/xin/GitHub_Weighted_lasso"$time
mkdir $backup_dir
cp script*.* $backup_dir

#cd $backup_dir
#git commit -m "update scripts"
#git push origin master
#cd $current_dir
