#!/bin/bash -l

err_report() {
  echo "Error on line $1 - script create_prefilt_Reads_subdir.sh"
  exit 2
}

trap 'err_report $LINENO' ERR

start=`date +%s.%N`
#---arguments---
wd=$(pwd)
fract=$1
reads_ext=$2
dataset=$3
Dts=$4
#maxjobs=$(echo $(nproc --all)/2.2 | bc)
maxjobs=$(echo $5/2.2 | bc)
#--- Main ----

mkdir -p $dataset/Reads/fraction_$fract

subsample_reads() {
  r=$1
  wd=$2
  dataset=$3
  fract=$4

  read_file=$(basename $r)
   if ! test -s $wd/$dataset/Reads/fraction_$fract/$read_file  #If file doesn't exit or if it exist but it is empty
      then
           all_lines_in=$(gzip -cd $r | wc -l)
           all_number_reads=$( echo $all_lines_in/4 | bc )
           subsample=$( echo $all_number_reads*$fract | bc )
           seqtk sample -s100 $r $subsample > $wd/$dataset/Reads/fraction_$fract/temponame_$read_file
           lines_in=$(cat $wd/$dataset/Reads/fraction_$fract/temponame_$read_file | wc -l)
           number_reads=$( echo $lines_in/4 | bc )
           gzip -c $wd/$dataset/Reads/fraction_$fract/temponame_$read_file > $wd/$dataset/Reads/fraction_$fract/$read_file
           rm $wd/$dataset/Reads/fraction_$fract/temponame_$read_file
           echo "      Subset $read_file created - Number of reads in subset: $number_reads"
   else
        echo "      Subset $read_file already created"
   fi

}

Rds=($(ls $wd/RAW_DATA/Reads/$Dts/*$reads_ext))

for r in "${Rds[@]}"
   do
      subsample_reads "$r" "$wd" "$dataset" "$fract" &
      if [[ $(jobs -r -p | wc -l) -ge $maxjobs ]]; then wait -n; fi
   done
wait

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
echo "INFO: Subsampling reads done in $runtime (s)"
