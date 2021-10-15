#!/bin/bash

skip_chunks=False
drop_analysis=False

#Command to run the pipeline on a server:
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -i|--imputed_target)
      IMPUTED="$2"
      shift # past argument
      shift # past value
      ;;
    -w|--wgs_target)
      WGS="$2"
      shift # past argument
      shift # past value
      ;;
    -bw|--bwgs_target)
      BWGS="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--threads)
      THREADS="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--skip_chunks)
      skip_chunks="True"
      shift # past argument
      shift # past value
      ;;
    -h|--help) 
      echo "Usage: imputation_accuracy.sh -i imputed_file.vcf.gz -w wgs_file.vcf.gz -bw wgs_HG00479.bcf.gz -t 10"
      exit 1 
      ;;
    -d|--drop_analysis)
      drop_analysis="True"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      echo "Unknown parameter passed: $1"
      echo "Usage: imputation_accuracy.sh -i imputed_file.vcf.gz -w wgs_file.vcf.gz"
      exit 1 
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

echo "Imputed VCF   : ${IMPUTED}"
echo "WGS VCF       : ${WGS}"
echo "WGS BCF       : ${BWGS}"
echo "Threads       : ${THREADS}"
echo "Skip chunks   : ${skip_chunks}"
echo "Skip analysis : ${drop_analysis}"

function max_bg_procs 
{
  if [[ $# -eq 0 ]] ; then
    echo "Usage: max_bg_procs NUM_PROCS.  Will wait until the number of background (&)"
    echo "           bash processes (as determined by 'jobs -pr') falls below NUM_PROCS"
    return
  fi
  local max_number=$((0 + ${1:-0}))
  while true; do
    local current_number=$(jobs -pr | wc -l)
    if [[ $current_number -lt $max_number ]]; then
      break
    fi
    sleep 1
  done
}


if [ $skip_chunks == "False" ]; then 
  #grab the header
  zcat ${IMPUTED} | head -10000 | grep '^#' > header
  #grab the non header lines
  zgrep -v '^#' ${IMPUTED} > variants
  #split into chunks with 1000 lines
  split -l 500000 variants
  #reattach the header to each and clean up
  for i in x*;do 
    cat header $i | bgzip > $i.vcf.gz
    tabix -p vcf $i.vcf.gz
    rm -f $i;
  done
  echo "Splitted Imputed file in chuncks of [500k]"
  rm -f header variants
  for i in x*.vcf.gz; do 
    filename=$(basename -- "$i")
    out_name="${filename%%.*}"
    bcftools view $i -Ob -o ${out_name}.bcf.gz --threads 4
    bcftools index ${out_name}.bcf.gz;
  done
  echo "BCF Imputed files Created"
fi

if [ $drop_analysis == "False" ]; then 
  n_file=`ls x*.vcf.gz | wc -l`
  #setup_scroll_area
  for splitted in x*.vcf.gz; do
      filename=$(basename -- "$splitted")
      file_name="${filename%%.*}"
      max_bg_procs $THREADS
      python3 /home/ec2-user/adriano/git/imputation_accuracy_calculator/Simpy.py --wgs ${WGS} --imputed ${file_name}.vcf.gz --bwgs ${BWGS} --bimputed ${file_name}.bcf.gz
  done


  echo "Joining files..."
  filename=$(basename -- "${IMPUTED}")
  out_name="${filename%%.*}"
  awk 'FNR>1 || NR==1' *_per_sample_results* > ${out_name}_per_sample_results.txt
  awk 'FNR>1 || NR==1' *_per_variant_results* > ${out_name}_per_variant_results.txt
  #delete tmp files
  echo "Deleting tmp files..."
  rm x*.vcf.gz*
  rm x*.bcf.gz*
  rm x*.vcf_per_variant_results*
  rm x*.vcf_per_sample_results*
fi

python3 /home/ec2-user/adriano/git/imputation_accuracy_calculator/rebuild_metrics.py --s ${out_name}_per_sample_results.txt

awk 'FNR>1 || NR==1' *_ImputationAccuracy.txt  > ${out_name}_per_sample_results.txt
rm *_ImputationAccuracy.txt
echo "Output-1    : ${out_name}_per_sample_results.txt"
echo "Output-2    : ${out_name}_per_variant_results.txt"
mv ${out_name}_per_sample_results.txt tmp_results
mv ${out_name}_per_variant_results.txt tmp_results