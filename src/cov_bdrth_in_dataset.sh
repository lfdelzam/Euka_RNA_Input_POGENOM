#!/bin/bash -l

err_report() {
    echo "Error on line $1 - script cov_bdrth_in_dataset.sh"
    exit 1
}

trap 'err_report $LINENO' ERR

# ---- arguments
mpileupfile=$1
bamfile=$2
outbamfile=$3
mag=$4
mincov=$5
minbreadth=$6
threads=$7
dataset=$8
samplename=$9
pdir="${10}"
m_ext="${11}"

wkd=$(pwd)

#--- Median coverage
cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

#---size
direct=$(echo "$dataset" | sed s/"_prefilt"//)

positions=$(awk 'BEGIN{i=0}; !/^>/ {i=i+length($0)} END {print i}' $wkd/RAW_DATA/Genomes/$direct/$mag$m_ext )

mkdir -p Genome_sizes
echo "genome size:" $positions > Genome_sizes/$mag.size

#---breadth
breadth=$(cut -f4 $mpileupfile | grep -cvw "0")
#breadthP=$(echo $breadth*100/$positions | bc -l )
breadthP=$( printf "%.2f \n" $(echo $breadth*100/$positions | bc -l ) )

if [ -f $wkd/RAW_DATA/gff_files/$direct/$mag.gff ]; then
  ExonBS=$( grep "exon" $wkd/RAW_DATA/gff_files/$direct/$mag.gff | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$5-$4+1}; END {print L}')
  echo "Exon bs:" $ExonBS >> Genome_sizes/$mag.size
  breadthExP=$( printf "%.2f \n" $(echo $breadth*100/$ExonBS | bc -l ) )
  echo "Genome: $mag - Sample: $samplename Median_coverage: $cov breadth (bs): $breadth - (% of genome): $breadthP - (Exons, %): $breadthExP"
else
  echo "Genome: $mag - Sample: $samplename Median_coverage: $cov breadth (bs): $breadth - (% of genome) $breadthP"
fi

mkdir -p 04_mergeable/$dataset/$pdir/$mag
bs_higherthan_min_cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk -v min=$mincov ' {if ($1 >= min) ++n} END { print n }')

#---selection of BAM files
if (( $(echo "$bs_higherthan_min_cov >= $minbreadth" | bc -l) )); then
 breadthratio=$( printf "%.2f \n" $(echo $bs_higherthan_min_cov*100/$breadth | bc -l ) )
 echo "Genome: $mag - Sample: $samplename - Number of bs with coverage higher than  $mincov : $bs_higherthan_min_cov - (% of breadth) $breadthratio"
 samtools sort -o $outbamfile --threads $threads $bamfile

fi
