#!/bin/bash 
#===================================
#  Bourne-Again shell script
#
#  Description:
#
#
#  Usage:
#
#
#  History:
#
#  Author:
#    Yi-Hsuan Chen
#    yihsuan@umich.edu
#===================================

#********************************************************
# Description:
#
#   Merge multiple images into a single image
#
#********************************************************

# data directory
wrkdir="./"

# output file name start
outfile_start=""

# dpi
dpi=300

# convert command
convert_command="convert"

usage="Usage: ./convert-merge_images.sh [-h|-v] [figures_to_merge] [figure_new]"

##################
# program start
##################

# read input parameters from command line
pram_num=$#

if [ $pram_num -eq 0 ]; then
  echo "No given file!"
  echo "program stop"
  exit 3
else

  for (( i=0; i<=$pram_num; i=i+1  ))
  do
    infile[i]=$1
    shift
  done
fi

# set append option
if [ ${infile[0]} == "-h" ]; then
  convert_command="convert +append"
elif [ ${infile[0]} == "-v" ]; then
  convert_command="convert -append"
else
  echo "ERROR: the first argument [${infile[0]}] must be [-h] (horizontal) or [-v] (vertical)"
  echo $usage
  exit
fi

# set the new file name
num_infile=$((${#infile[@]}-2))
new_file=${infile[$num_infile]}
files=""

#echo $num_infile
#echo $convert_command
#echo $new_file

# check with user
echo "------------------------"
echo "Merge multiple images into a single image"
echo ""
echo "Process files"
echo ""
for ((i=1; i<$num_infile; i=i+1))
do
  nn1=$((num_infile-1))
  file1=${infile[$i]}
  if [ ! -f $file1 ]; then
    echo "ERROR: file [$file1] does not exist"
    exit 1
  fi
  files="$files $file1"
  echo "  input files        : $i/$nn1, [${file1}]"
done
  echo ""
  #echo "files: [$files]"
  echo "  output flie        : [$new_file]"
  echo ""
  echo "  convert command    : [$convert_command]"

echo "------------------------"
read -p "Is is correct? (y/n)  " choice
echo " "
if [ ! $choice ] || [ ! $choice == "y" ]; then
  echo "Cancel by user"
  echo "program stop"
  exit 5
fi

# use convert
$convert_command $files $new_file && \
   echo "Done. create [$new_file]" || exit 5

exit 0
