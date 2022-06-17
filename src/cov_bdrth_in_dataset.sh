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

wkd=$(pwd)

#--- Median coverage
cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

#---size
check=$(echo "$dataset" | cut -d "_" -f2)

if [[ "$check" == "prefilt" ]]; then
   direct=$(echo $dataset | sed s/"_prefilt"//)
else
    direct="$dataset"
fi

positions=$( grep -v ">" $wkd/RAW_DATA/Genomes/$direct/$mag.* | tr -d "\n" | wc -c)

mkdir -p Genome_sizes
echo "genome size:" $positions > Genome_sizes/$mag.size

#---breadth
#non_zero=$(cut -f4 $mpileupfile | grep -cvw "0")
#breadth=$(echo $non_zero*100/$positions | bc -l )
breadth=$(cut -f4 $mpileupfile | grep -cvw "0")

if [ -f $wkd/RAW_DATA/gff_files/$direct/$mag.gff ]; then
  ExonBS=$( grep "exon" $wkd/RAW_DATA/gff_files/$direct/$mag.gff | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$5-$4+1}; END {print L}')
  echo "Exon bs:" $ExonBS >> Genome_sizes/$mag.size
  breadthExP=$(echo $breadth*100/$ExonBS | bc -l )
  echo "Genome:" $mag "- Sample:" $samplename "Median_coverage:" $cov " breadth (bs):" $breadth "breadth (Exons, %):" $breadthExP
else
  echo "Genome:" $mag "- Sample:" $samplename "Median_coverage:" $cov " breadth (bs):" $breadth
fi

mkdir -p 04_mergeable/$dataset/$pdir/$mag
bs_higherthan_min_cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk -v min=$mincov ' {if ($1 >= min) ++n} END { print n }')

#---selection of BAM files and subsample -old
#if (( $(echo "$breadth >= $minbreadth" | bc -l) )) && (( $(echo "$cov >= $mincov" | bc -l) )); then
#  limite=$(echo "scale=3; $mincov/$cov" | bc )
#  samp=$(echo "scale=3; ($limite)+10" | bc)
#  samtools view -Sbh --threads $threads -s $samp $bamfile | samtools sort -o $outbamfile --threads $threads

#---selection of BAM files
if (( $(echo "$bs_higherthan_min_cov >= $minbreadth" | bc -l) )); then

 samtools sort -o $outbamfile --threads $threads $bamfile

fi
