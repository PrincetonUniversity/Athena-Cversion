#!/bin/bash

# Script to join ppm files
# Last updated August 24, 2012 by Sasha Philippov

# Usage is, e.g.,
# joinppm.sh -i FP -o comb/FPc -f 1:4 -p 16 -s d.ppm -x 2 -y 2

# That will take the d.ppm files from the id*/ sub-directories of the current
# directory and put the output (changing the filename from FP*.ppm to FPc*.ppm
# so you don't confuse the two) in the sub-directory comb (which it will
# create, if necessary). -x and -y specify the mpi decomposition of the
#computational domain. See also below.

scrdir=$(dirname $0)
script=$(basename $0)

# set defaults
dir=$(pwd)
ext="ppm"
existid0=1
usesub=1
notify=0
eight=0
multi=0
waitfor=0
sdt=120
fs=1
ret=0

# ***************************************************************************

usage ()
{
  echo "This script joins ppm files."
  echo "      (e-mail lemaster@astro for assistance)"
  echo ""
  echo "Options:"
  echo "  -i inbase"
  echo "          Base of input files (e.g. Turb for Turb-id1.0000.ppm)"
  echo "  -o outbase"
  echo "          Base of output files (e.g. combined for combined.0000.ppm)"
  echo "  [-s extension]"
  echo "          File extension (e.g. d.ppm for Turb.0000.d.ppm) [ppm]"
  echo "  -f {f|f0:fn|f0:fn:fi}"
  echo "          Join specific file number (f) OR file range (f0 thru fn) OR"
  echo "          sequence of files (f0 thru fn at intervals fi)"
  echo "  -p {nproc|p0:pn}"
  echo "          Number of processors (nproc) OR processor range (p0 thru pn)"
  echo "  [-d dir]"
  echo "          Input directory [.]"
  echo "  [-n]"
  echo "          E-mail completion notification to $(whoami)@astro.princeton.edu"
  echo "  [-e]"
  echo "          Join into eight ppm slabs (a thru h) in z-dir instead of one"
  echo "  [-w]"
  echo "          Wait while files are transfered by another script"
  echo "  [-z]"
  echo "          Root process does not write to id0 sub-directory"
echo "  -x Nx"
echo "          Number of zones in x direction (Nx)"
echo "  -y Ny"
echo "          Number of zones in y direction (Ny)"
}

# ***************************************************************************

checkparset ()
{
  if [ "X$2" == "X" ]; then
    echo "Error:  $1"
    usage
    exit 1
  fi
}

# ***************************************************************************

checkopts ()
{ # validate parameters
  checkparset "Input base name not set! (use -i)" ${ibase}
  checkparset "Output base name not set! (use -o)" ${obase}

  checkparset "Processor range not set! (use -p)" ${p0}
  checkparset "Processor range not set! (use -p)" ${pn}

  if [ ${p0} -lt 0 ]; then
    echo "The starting processor number cannot be negative."
    usage
    exit 1
  fi
  if [ ${pn} -lt 0 ]; then
    echo "The number of processors must be positive."
    usage
    exit 1
  fi
  if [ ${pn} -lt ${p0} ]; then
    echo "The ending processor number cannot be smaller than the starting."
    usage
    exit 1
  fi
  if [ ${eight} -ne 0 ]; then
    nproc=$[${pn}+1]
    if [ ${p0} -ne 0 ]; then
      echo "Combining into eight files only works when starting with proc 0."
      usage
      exit 1
    fi
    rem=$[${nproc} % 8]
    if [ ${rem} -ne 0 ]; then
      echo "The number of processors must be divisible by eight."
      usage
      exit 1
    fi
    div=$[${nproc} / 8]
  else
    if [ $[${pn}-${p0}+1] -gt 512 ]; then
      multi=1
      div=512
    fi
  fi

  checkparset "File range not set! (use -f)" ${f0}
  checkparset "File range not set! (use -f)" ${fn}
  checkparset "File range not set! (use -f)" ${fs}

  if [ ${f0} -lt 0 ]; then
    echo "The starting file number cannot be negative."
    usage
    exit 1
  fi
  if [ ${fn} -lt ${f0} ]; then
    echo "The ending file number cannot be smaller than the starting."
    usage
    exit 1
  fi
  if [ ${fs} -lt 0 ]; then
    echo "The file interval cannot be negative."
    usage
    exit 1
  fi
}

# ***************************************************************************

padstring ()
{
  # pad the file number with zeros
  if [ $FNO -lt 10 ]; then
    FNOSTR="000${FNO}"
  elif [ $FNO -lt 100 ]; then
    FNOSTR="00${FNO}"
  elif [ $FNO -lt 1000 ]; then
    FNOSTR="0${FNO}"
  else
    FNOSTR="${FNO}"
  fi
}

sourcelist ()
{
  # construct list of source files
  xp0=$1
  xpn=$2
  xibase=$3
  xcap=$4
  xusesub=$5
  xext=$6

  xp1=${xp0}
  src=""
  sub=""
  if [ ${existid0} -eq 0 ]; then
    if [ ${xp0} -eq 0 ]; then
      src="${xibase}.${FNOSTR}.${xext}"
      xp1=1
    fi
  fi
  for ((a=${xp1} ; a<=${xpn} ; a++)) ; do
    if [ ${xusesub} -ne 0 ]; then
      sub="id${a}/"
    fi
    if [ ${a} -ne 0 ]; then
      src="${src} ${sub}${xibase}${xcap}${a}.${FNOSTR}.${xext}"
    else
      src="${src} ${sub}${xibase}.${FNOSTR}.${xext}"
    fi
  done
}

joinlist ()
{
  outfile="$6.${FNOSTR}.$7"
  echo "Creating ${outfile}..."

  sourcelist $1 $2 $3 $4 $5 $7

  # join the files
montage ${src} -tile ${Nx}x${Ny} -geometry \>0+0 -flip ${outfile}
convert ${outfile} -flip ${outfile}
    xret=$?
  if [ ${xret} -ne 0 ]; then
    rm ${outfile} 2>/dev/null
    echo "Error joining ppm file ${outfile}"
    if [ ${notify} -ne 0 ]; then
      echo "${dir}/${outfile}" | mail $(whoami)@astro.princeton.edu -s "Error joining ppm file ${outfile}"
    fi
    exit ${xret}
  fi
  return ${xret}
}

# ***************************************************************************

while getopts "enwzd:f:b:i:o:p:s:x:y" opt ; do
  case $opt in
    "e" )
#      echo "writing eight ppms a-h"
      eight=1
      ;;
    "n" )
#      echo "notify of completion"
      notify=1
      ;;
    "w" )
#      echo "will wait for source files"
      waitfor=1
      ;;
    "z" )
#      echo "root process doesn't write to id0"
      existid0=0
      ;;
    "d" )
      dir=${OPTARG}
#      echo "source file full path ${dir}"
      ;;
    "f" )
      cnt=$(echo "${OPTARG}" | awk "{ cnt = split(\$0,a,\":\"); print cnt }")
      if [ ${cnt} -eq 1 ]; then
        f0=${OPTARG}
        fn=${f0}
#        echo "file ${f0}"
      else
        if [ ${cnt} -eq 2 ]; then
          f0=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[1] }")
          fn=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[2] }")
#          echo "files ${f0} through ${fn}"
        else
          f0=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[1] }")
          fn=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[2] }")
          fs=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[3] }")
#          echo "files ${f0} through ${fn} with interval ${fs}"
	fi
      fi
      ;;
    "b" )
      ibase=${OPTARG}
#      echo "source file base ${ibase}"
      ;;
    "i" )
      ibase=${OPTARG}
#      echo "source file base ${ibase}"
      ;;
    "o" )
      obase=${OPTARG}
#      echo "output file base ${obase}"
      ;;
    "p" )
      cnt=$(echo "${OPTARG}" | awk "{ cnt = split(\$0,a,\":\"); print cnt }")
      if [ ${cnt} -eq 1 ]; then
        p0=0
        pn=$[${OPTARG}-1]
#        echo "procs 0 through ${pn}"
      else
        p0=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[1] }")
        pn=$(echo "${OPTARG}" | awk "{ split(\$0,a,\":\"); print a[2] }")
#        echo "procs ${p0} through ${pn}"
      fi
      ;;
    "s" )
      # e.g. ppm or d.ppm
      ext=${OPTARG}
#      echo "source file extension ${ext}"
      ;;
    "x" )
      Nx=${OPTARG}
#      echo "output file base ${obase}"
      ;;

    "y" )
      Ny=${OPTARG}
#      echo "output file base ${obase}"
    ;;

    * )
      badsyn=1
#      echo "invalid option"
  esac
done
shift $(($OPTIND - 1))

# validate parameters
checkopts

# change to the source directory
if [ "X${dir}" != "X" ]; then
  if [ ${waitfor} -ne 0 ]; then
    until [ -d $4 ]; do
      sleep ${sdt}
    done
  fi
  cd ${dir}/
  if [ $? -ne 0 ]; then
    echo "Error: Couldn't cd to ${dir}."
    exit 1
  fi
  pwd
fi

# wait for the source files to be copied by some other script
if [ ${waitfor} -ne 0 ]; then
  FNO=${fn}
  padstring

  if [ ${usesub} -ne 0 ]; then
    wfile="id${pn}/${ibase}-id${pn}.${FNOSTR}.${ext}"
  else  
    wfile="${ibase}-id${pn}.${FNOSTR}.${ext}"
  fi
  echo "Waiting for ${wfile} to be created by external script..."

  cnt=$[$(ls -1 ${wfile} | wc -l 2>/dev/null)]
  while [ ${cnt} -eq 0 ]; do
    sleep ${sdt}
    cnt=$[$(ls -1 ${wfile} | wc -l 2>/dev/null)]
  done
fi

# create the output directory if needed
osub=$(dirname ${obase})
if [ "X${osub}" != "X" ]; then
  mkdir ${osub}/ 2>/dev/null
fi

# for FNO in `seq ${f0} ${fn}`; do
for ((FNO=${f0} ; FNO <= ${fn} ; FNO += ${fs})) ; do
  # pad the file number with zeros
  padstring

  # construct list of source files and join them
  if [ ${eight} -eq 0 ]; then
    # write to single ppm file
    if [ ${multi} -eq 0 ]; then
      joinlist ${p0} ${pn} ${ibase} "-id" ${usesub} ${obase} ${ext}
      ret=$[${ret}+$?]
    else
      # too many files to join all at once; do in stages instead
      ret1=0
      cnt=0
      pna=$[$[$cnt+1]*${div}-1]
      while [ ${pna} -lt ${pn} ]; do
        p0a=$[${cnt}*${div}]
        pna=$[$[$cnt+1]*${div}-1]
        if [ ${pna} -gt ${pn} ]; then
          pna=${pn}
        fi
        cnt=$[${cnt}+1]
        joinlist ${p0a} ${pna} ${ibase} "-id" ${usesub} ${obase}-${cnt} ${ext}
        ret1=$[${ret1}+$?]
      done
      if [ ${ret1} -eq 0 ]; then
        joinlist 1 ${cnt} ${obase} "-" 0 ${obase} ${ext}
        ret1=$?
        if [ ${ret1} -eq 0 ]; then
          for p0a in `seq 1 ${cnt}`; do
            rm ${obase}-${p0a}.${FNOSTR}.${ext}
          done
        fi
      fi
      ret=$[${ret}+${ret1}]
    fi
  else
    # write to eight ppm slabs (a thru h)
    joinlist $[0*${div}] $[1*${div}-1] ${ibase} "-id" ${usesub} ${obase}-a ${ext}
    ret=$[${ret}+$?]
    joinlist $[1*${div}] $[2*${div}-1] ${ibase} "-id" ${usesub} ${obase}-b ${ext}
    ret=$[${ret}+$?]
    joinlist $[2*${div}] $[3*${div}-1] ${ibase} "-id" ${usesub} ${obase}-c ${ext}
    ret=$[${ret}+$?]
    joinlist $[3*${div}] $[4*${div}-1] ${ibase} "-id" ${usesub} ${obase}-d ${ext}
    ret=$[${ret}+$?]
    joinlist $[4*${div}] $[5*${div}-1] ${ibase} "-id" ${usesub} ${obase}-e ${ext}
    ret=$[${ret}+$?]
    joinlist $[5*${div}] $[6*${div}-1] ${ibase} "-id" ${usesub} ${obase}-f ${ext}
    ret=$[${ret}+$?]
    joinlist $[6*${div}] $[7*${div}-1] ${ibase} "-id" ${usesub} ${obase}-g ${ext}
    ret=$[${ret}+$?]
    joinlist $[7*${div}] $[8*${div}-1] ${ibase} "-id" ${usesub} ${obase}-h ${ext}
    ret=$[${ret}+$?]
  fi
done

if [ ${notify} -ne 0 ]; then
  if [ ${ret} -eq 0 ]; then
    echo "${dir}" | mail $(whoami)@astro.princeton.edu -s "Done joining, no errors: ${dir}"
  else
    echo "${dir}" | mail $(whoami)@astro.princeton.edu -s "Done joining, with errors: ${dir}"
  fi
fi

exit ${ret}
