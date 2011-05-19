#!/bin/bash
minuso=$1
outbasename=$2
inbasename=$3
s=$4
e=$5

homedir="/aether1/askinner/svnathena_r1317/vis/matlab/join_bin"

usage_err="Usage:  join_bin.sh -o <OUTBASENAME> <INBASENAME> <STARTCYCLE> <ENDCYCLE>"

if [[ $minuso != "-o" ]]
  then echo $usage_err 
  exit
elif (( "$s" < "0" ))
  then echo $usage_err 
  exit
elif (( "$e" < "$s" ))
  then echo $usage_err 
  exit
fi

gcc -o $homedir/join_bin $homedir/join_bin.c -lm

nproc=`ls -ld id* | wc -l`
nproc=$((nproc-1))
echo $nproc
for c in `seq $s $e`
do
  num=$(printf "%04d" "$c")
  arglist="id0/$inbasename.$num.bin"

  for k in `seq 1 $nproc`
  do
    arglist="$arglist id$k/$inbasename-id$k.$num.bin"
  done

  outfilename=$outbasename.$num.bin
  echo "Cycle $c:"
  echo "  Stitching $((nproc+1)) processors together to form $outfilename"
#  echo $arglist
  $homedir/join_bin -o $outfilename $arglist
done



