#!/bin/bash -l
wd=$(pwd)
err_report() {
    echo "Error on line $1 - script Euka_RNA_Input_POGENOM.sh"
    cd $wd
    if [ -z "$mag" ]; then mag=$(fullname=$(basename $(ls RAW_DATA/Genomes/$dataset/*$genomes_ext | head -1)) ; echo ${fullname%.*}); fi
    if [ -z "$samples" ]; then samples=$(fulln=$(basename $(ls RAW_DATA/Reads/$dataset/*$fwd_index$reads_ext | head -1)) ; echo ${fulln%$fwd_index$reads_ext}); fi
    mess1="Check if the 'RNA_Euka_ip_env' has been activated - command:\n    'conda activate RNA_Euka_ip_env'"
    mess2="Use the command:\n    'snakemake -s snakefiles/Euka_RNA_step1_pogenom_input --unlock'"
    mess3="Use the command:\n    'snakemake -s snakefiles/Euka_RNA_step_pogenom_input --config my_mag='"$mag"' my_samples='"$samples"' --unlock'"
    mess="or\n    'snakemake -s snakefiles/Euka_RNA_step2 --config mag_name='"$mag"' --unlock' or\n    'snakemake -s snakefiles/Euka_RNA_step2_B --config mag_name='"$mag"' --unlock'"
    if [ "$1" == 77 ]; then echo -e "if you are using conda, check if the 'RNA_Euka_ip_env' has been activated - command:\n    'conda activate RNA_Euka_ip_env' "; fi
    if [ "$1" == 79 ]; then echo -e "TIP 1 - Look at $wd/log_files/samples_filter_$dataset.log\nTIP 2 - $mess1\nTIP 3 - Use the command:\n    'snakemake -s snakefiles/Euka_RNA_step_filter --unlock'\n     and run the pipeline again"; fi
    if [ "$1" == 95 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag.coverage_breadth.log\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 98 ] || [ "$1" == 101 ] ; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag.prefilt_readcounts.log\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 108 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset.$mag"_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess3 $mess\n     and run the pipeline again"; fi
    if [ "$1" == 124 ]; then echo -e "TIP 1 - Look at\n    $wd/log_files/$dataset.$mag.coverage_breadth.log or\n    $wd/log_files/$dataset.$mag"_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess3\n     and run the pipeline again"; fi
    if [ "$1" == 129 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_coverage_breadth.log"\nTIP 2 - $mess1\nTIP 3 - $mess2\n     and run the pipeline again"; fi
    if [ "$1" == 131 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_vcf_files.log"\nTIP 2 - $mess1\nTIP 3 - $mess2 $mess\n     and run the pipeline again"; fi
    if [ "$1" == 140 ] || [ "$1" == 143 ]; then echo -e "TIP 1 - Look at $wd/log_files/$dataset"_Genomes_readcounts.log"\nTIP 2 - $mess1\nTIP 3 - $mess2\n     and run the pipeline again\nTIP 4 - If you selected read_counts_per_gene : yes, check if all the genomes in the dataset have its gff file"; fi
    if test -f "temporal"; then rm temporal; fi
    exit 1
}
trap 'err_report $LINENO' ERR
#Default options
configFile=$wd/config_files/Euka_RNA_Input_POGENOM_config.json
#----Argument parse----------------------
for a in "$@"
do
case $a in
  -d=*|--path_to_config_file=*)
  if [ -z "${a#*=}" ];  then
  echo "value to argument -d No supplied"
  exit 0
  else configFile="${a#*=}"
  fi
  shift # past argument
  ;;

  *)
  echo -e "\nUsage: bash Euka_RNA_Input_POGENOM.sh [options]\n -d=<absolute path to configFile. Default=$configFile>\n"
  echo -e 'Description:\nThis program executes a pipeline that generates the required input files for POGENOM.\nThe aim of this pipeline is to increase the reproducibility of the data analysis, and to simplify the use of POGENOM.\nPOGENOM is a computer program that calculates several population genetic parameters for a genome in relation to a set of samples (https://github.com/EnvGen/POGENOM).'
  exit 0
  ;;

esac
done
if [[ "$configFile" != /* ]] || [ -z "$configFile" ]; then
    echo "Please provide an absoltute path to configfile e.g., bash Euka_RNA_Input_POGENOM.sh '/absolute/path/to/configfile' "
    exit 0
fi

cat $configFile | sed s/"[{|}]"//g | sed s/":"/"="/g | sed s/",$"//g | sed s/" ="/"="/g | sed s/"= "/"="/g | sed s/'"'//g | sed s/" "//g > temporal
. temporal

if [[ "$workdir" != /* ]] || [ -z "$workdir" ]; then
    echo "Please provide an absoltute path to the working directory in configfile e.g., 'workdir': '/absolute/path/to/working_directory/' "
    echo "Please double-check the absolute path to working directory $workdir and to configfile $configFile"
    rm temporal
    exit 0
fi
#Checking key parameters setting
options=("$dataset" "$min_coverage" "$min_breadth" "$min_bsq_for_cov_median_calculation" "$threads" "$genomes_ext" "$reads_ext" "$fwd_index" "$rev_index" "$hisat2_params" "$htseq_params" "$mapqual" "$freebayes_parameters" "$vcffilter_qual")
for o in "${options[@]}"; do if [ -z "$o" ]; then echo "A key parameter is undefined, please check in the config_files/Euka_RNA_Input_POGENOM_config.json file the parameters used"; exit 1; fi; done
if [[ $snakemake_extra_params == *","* ]]; then extra_params=$( echo $snakemake_extra_params | sed s/","/" "/g); else extra_params=$snakemake_extra_params; fi
echo "INFO: Starting Euka_RNA_Input_POGENOM pipeline - Working directory: $workdir"
#----Using prefilt mode - full workflow
mkdir -p $workdir/log_files
if  [[ "$mode_prefilt" == TRUE ]]; then
#Checking key parameters setting
options2=("$fraction" "$temp_sub_Reads_dir")
for p in "${options2[@]}"; do if [ -z "$p" ]; then echo 'A key parameter in "mode_prefilt" is undefined, please check in the config_files/Euka_RNA_Input_POGENOM_config.json file the parameters used'; exit 1; fi; done
if [[ "$read_counts_per_gene" == TRUE ]] && [ -z "$htseq_params" ]; then echo 'The key parameter "htseq_params" is undefined, please check in the config_files/Euka_RNA_Input_POGENOM_config.json file the parameters used'; exit 1; fi
# main - mode prefilt
         cd $workdir
         echo "INFO: Generating Reads subsets - Fraction used $fraction"
         bash src/create_prefilt_Reads_subdir.sh $fraction $genomes_ext $reads_ext $temp_sub_Reads_dir $dataset
         echo "INFO: Calculating Genome coverage - sub-samples - Coverage threshold $min_coverage"
         snakemake -s snakefiles/Euka_RNA_step_filter -j $threads $extra_params 2>log_files/samples_filter_$dataset.log
         if [[ "$remove_subreads" == TRUE ]] && test -d "$temp_sub_Reads_dir/Reads"; then
              echo "WARNING: You have chosen to remove $temp_sub_Reads_dir/Reads/"
              rm -rf $temp_sub_Reads_dir/Reads/
         fi
         result_dir="PREFILT/"$dataset"/params_cov_"$min_coverage"_mpq_"$mapqual"_bq_"$min_bsq_for_cov_median_calculation"_fr_"$fraction
        # if [ -s $result_dir/Selected_samples_Genomes.txt ]; then
             file_empty=$(grep -v "#" $result_dir/Selected_samples_Genomes.txt | wc -l)
             if [ "$file_empty" -eq 0 ]; then
                echo -e "INFO: With the current parameter setting: Dataset $dataset - Fraction $fraction - Coverage threshold $min_coverage - Min-base quality $min_bsq_for_cov_median_calculation - Mapping quality $mapqual\n      There is no Genome - sample with Estimated Median Coverage higher than threshold.\n      A vcf file cannot be created\n"
             else
                grep -v "#" $result_dir/Selected_samples_Genomes.txt | while read line
                do
                  mag=$(echo $line | cut -d " " -f1)
                  samples=$(echo $line | cut -d " " -f2)
<<<<<<< HEAD
                  echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Genome: $mag - Coverage threshold: $min_coverage - Breadth threshold: $min_breadth "
=======
                  echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Genome: $mag - Median coverage threshold: $min_coverage - Breadth threshold: $min_breadth "
>>>>>>> a8cad084431f54cb8fcc2e6609cb7bb0b0fd4541
                  snakemake -s snakefiles/Euka_RNA_step_pogenom_input step1_all --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag.coverage_breadth.log
                  if [[ "$read_counts_per_gene" == TRUE ]] && [ -f $workdir/RAW_DATA/gff_files/$dataset/$mag.gff ]; then
                     if [ ! -z "$( ls -A 04_mergeable/"$dataset"_prefilt/params_cov_"$min_coverage"_bdth_"$min_breadth"_mpq_"$mapqual"_bq_"$min_bsq_for_cov_median_calculation"/$mag/*.bam)" ]; then
                         snakemake -s snakefiles/Euka_RNA_step_pogenom_input readcounts_mergeable_all --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag.prefilt_readcounts.log
                         echo "INFO: Counting reads analysis has been performed using the BAM file(s) in 04_mergeable/"$dataset"_prefilt/params_cov_"$min_coverage"_bdth_"$min_breadth"_mpq_"$mapqual"_bq_"$min_bsq_for_cov_median_calculation"/$mag/"
                     else
                         snakemake -s snakefiles/Euka_RNA_step_pogenom_input readcounts_all --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag.prefilt_readcounts.log
                         echo "INFO: Counting reads analysis has been performed using the all BAM files in 02_MAPPING/"$dataset"_prefilt/$mag/"
                     fi
                  else
                    echo "Warning: Counting reads analysis was no realised on Genome $mag"; if [ $read_counts_per_gene == TRUE ]; then echo "         please provide the required gff file: $workdir/RAW_DATA/gff_files/$dataset/$mag.gff"; fi
                  fi
                  echo "INFO: Generating VCF files - Genome $mag"
                  snakemake -s snakefiles/Euka_RNA_step_pogenom_input vcf --config my_mag="$mag" my_samples="$samples" -j $threads $extra_params 2> log_files/$dataset.$mag"_vcf_files.log"
                done
             fi
             no_genome=$(grep "#" $result_dir/Selected_samples_Genomes.txt | wc -l)
             if [ "$no_genome" -ne 0 ]; then
                 echo "**********************************************"
                 echo "The following Genome(s) has(have) not been analysed"
                 grep "#" $result_dir/Selected_samples_Genomes.txt
                 echo -e "**********************************************\n"
             fi
         echo "INFO: Euka_RNA_Input_POGENOM pipeline is done !!!"
#         else
#           echo -e "ERROR: file $result_dir/Selected_samples_Genomes.txt was not created\nTIP 1 - Look at $wd/log_files/samples_filter_$dataset.log\nTIP 2 - Check if the 'RNA_Euka_ip_env' has been activated - command:\n    'conda activate RNA_Euka_ip_env' "
#         fi
rm temporal
exit 0
fi
#---End of mode prefilt
#---Option when analysing a dataset without prefilt
cd $workdir
<<<<<<< HEAD
echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Coverage threshold: $min_coverage - Breadth threshold: $min_breadth "
=======
echo "INFO: Calculating Genome Median coverage and breadth - Dataset: $dataset - Median coverage threshold: $min_coverage - Breadth threshold: $min_breadth "
>>>>>>> a8cad084431f54cb8fcc2e6609cb7bb0b0fd4541
   snakemake -s snakefiles/Euka_RNA_step1_pogenom_input step1_all -j $threads $extra_params 2> log_files/$dataset"_Genomes_coverage_breadth.log"
echo "INFO: Generating VCF files"
   snakemake -s snakefiles/Euka_RNA_step1_pogenom_input vcf -j $threads $extra_params 2> log_files/$dataset"_Genomes_vcf_files.log"
if [[ "$read_counts_per_gene" == TRUE ]]; then
    mergeable_dir="$workdir/04_mergeable/$dataset/params_cov_"$min_coverage"_bdth_"$min_breadth"_mpq_"$mapqual"_bq_"$min_bsq_for_cov_median_calculation
    MAGs=( $(ls -d $mergeable_dir/* ) )
    for m in "${MAGs[@]}"
      do
        mag=$(basename $m)
        if [ -f $workdir/RAW_DATA/gff_files/$dataset/$mag.gff ]; then
           if [ ! -z "$( ls -A $mergeable_dir/$mag/*.bam)" ]; then
               snakemake -s snakefiles/step1_read_counts all_m --config mag_name="$mag" -j $threads $extra_params 2> log_files/$dataset.$mag.readcounts.log
               echo "INFO: Counting reads analysis has been performed using the BAM file(s) in $mergeable_dir/$mag/"
           else
               snakemake -s snakefiles/step1_read_counts all --config mag_name="$mag" -j $threads $extra_params 2> log_files/$dataset.$mag.readcounts.log
               echo "INFO: Counting reads analysis has been performed using the all BAM files in 02_MAPPING/"$dataset"_prefilt/$mag/"
           fi
        else
          echo "Warning: Counting reads analysis was no realised on Genome $mag"
          echo "         please provide the required gff file: $workdir/RAW_DATA/gff_files/$dataset/$mag.gff"
        fi
      done
fi

    rm temporal

echo 'INFO: Euka_RNA_Input_POGENOM pipeline is done !!!'
