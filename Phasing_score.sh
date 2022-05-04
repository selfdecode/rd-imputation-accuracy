#!/bin/bash
if [ -z "$*" ]; then 
    echo "
author: SELFDECODE
contact email : adriano@selfdecode.com

SelfDecode software to analyze multiple phasing and imputation softwares

Parameters:
"
      echo "      -h|--help
              show help"
      echo "      -i|--input
              input file"
      echo "      -r|--reference
              full path reference file without extension."
      echo "      -t|--threads
              number of cpus to use."
      echo "      -o|--output
              output prefix. No extension."
      echo "      -c|--chr
              chomosome to analyze allowed [1-22 X]. NO chr prefix. Remove also prefix from VCF"
      echo "      -pbeagle|--phase_beagle
              skip Beagle haplotype estimation"
      echo "      -shapeit|--shapeit
              skip ShapeIT Phasing"
      echo "      -eagle|--eagle
              skip Eagle phasing"
      echo "      -bigref|--BIGREF
              use this option if you get memory allocate error during accuracy evaluation"

echo "
[base] Usage: ./Phasing_score.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>
[skip] Usage: ./Phasing_score.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -shapeit no 
[memo] Usage: ./Phasing_score.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -bigref yes

"
    exit 1
fi

BEAGLE_PHASE=true
EAGLE=true
SHAPEIT=true
BIGREF=false


POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -i|--input)
      INPUT="$2"
      shift 
      shift
      ;;
    -r|--reference)
      REFERENCE="$2"
      shift 
      shift
      ;;
    -t|--threads)
      THREADS="$2"
      shift 
      shift
      ;;
    -c|--chr)
      CHROMOSOME="$2"
      shift 
      shift
      ;;
    -pbeagle|--phase_beagle)
      BEAGLE_PHASE="false"
      shift 
      shift
      ;;
    -eagle|--eagle)
      EAGLE="false"
      shift 
      shift
      ;;
    -shapeit|--shapeit)
      SHAPEIT="false"
      shift 
      shift
      ;;
    -bigref|--BIGREF)
      BIGREF="true"
      shift 
      shift
      ;;
    -o|--output)
      OUTPUT="$2"
      shift 
      shift
      ;;
    -h|--help) 
      echo "

author: SELFDECODE
contact email : adriano@selfdecode.com

SelfDecode software to analyze multiple phasing and imputation softwares

Parameters:
"
      echo "      -h|--help
              show help"
      echo "      -i|--input
              input file"
      echo "      -r|--reference
              full path reference file without extension."
      echo "      -t|--threads
              number of cpus to use."
      echo "      -o|--output
              output prefix. No extension."
      echo "      -c|--chr
              chomosome to analyze allowed [1-22 X]. NO chr prefix. No chr prefix in VCF file"
      echo "      -pbeagle|--phase_beagle
              skip Beagle haplotype estimation"
      echo "      -shapeit|--shapeit
              skip ShapeIT Phasing"
      echo "      -eagle|--eagle
              skip Eagle phasing"
      echo "      -bigref|--BIGREF
              use this option if you get memory allocate error during accuracy evaluation"

echo "
[base] Usage: ./Phasing_score.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>
[skip] Usage: ./Phasing_score.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -ibeagle no -impute5 no
[memo] Usage: ./Phasing_score.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -bigref yes

"
      exit 1 
      ;;
    *)    # unknown option
      echo "Unknown parameter passed: $1"
      exit 1 
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

echo "Input         : ${INPUT}"
echo "Output        : ${OUTPUT}"
echo "Reference     : ${REFERENCE}"
echo "Chromosome    : ${CHROMOSOME}"
echo "Threads       : ${THREADS}"
echo "BEAGLE PHASE  : ${BEAGLE_PHASE}"
echo "EAGLE         : ${EAGLE}"
echo "SHAPEIT       : ${SHAPEIT}"
echo "BIGREF        : ${BIGREF}"

#######################################################
#######################################################
#SOFTWARE PATH (EDITABLE PART)
time="/usr/bin/time -f"
eagle2="/home/ec2-user/adriano/imputation/phase2/software/eagle2.4.1/Eagle_v2.4.1/eagle"
shapeit4="/home/ec2-user/adriano/imputation/phase2/software/shapeit4/shapeit4-4.2.1/bin/shapeit4.2"
beagle5="java -jar /home/ec2-user/adriano/imputation/phase2/software/beagle5.2/beagle.25Mar22.4f6.jar"
bref3="java -Xmx8g -jar /home/ec2-user/adriano/imputation/phase2/software/beagle5.2/bref3.29May21.d6d.jar"
simpy="/home/ec2-user/adriano/git/rd-imputation-accuracy/bin/Simpy.py"
imputation_accuracy="/home/ec2-user/adriano/git/rd-imputation-accuracy/imputation_accuracy.sh"
base_output="/home/ec2-user/adriano/imputation/phase2/results"
#GENETIC RECOMBINATION MAP PATH
map_eagle2="/home/ec2-user/adriano/imputation/phase2/software/eagle2.4.1/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
map_beagle5="/home/ec2-user/adriano/imputation/phase2/genetic_map/plink.chr${CHROMOSOME}.GRCh38.map"
map_shapeit4="/home/ec2-user/adriano/imputation/phase2/genetic_map/chr${CHROMOSOME}.b38.gmap.gz"
#WGS PATH for accuracy
wgs_subset="/home/ec2-user/adriano/imputation/phase2/wgs/20.reference_panel.30x.hg38.540children_phased.vcf.gz"
bwgs_subset="/home/ec2-user/adriano/imputation/phase2/wgs/20.reference_panel.30x.hg38.540children_phased.vcf.gz"
#######################################################
#######################################################

if [ ! -e "${wgs_subset}" ]; then
  echo "CREATION OF WGS SUBSET IN VCF FORMAT"
  #echo "chr${CHROMOSOME} ${CHROMOSOME}" > chr${CHROMOSOME}_${CHROMOSOME}.txt
  #bcftools annotate --rename-chrs chr${CHROMOSOME}_${CHROMOSOME}.txt ${wgs_subset_chr} \
  #  -Oz -o ${wgs_subset} \
  #  --thread ${THREADS} && tabix ${output_name} -f
  bcftools view ${wgs_subset} -Ob -o ${bwgs_subset} && bcftools index ${bwgs_subset}
  # rm chr${CHROMOSOME}_${CHROMOSOME}.txt
fi

#DATA PATH
if [ ! -e "${INPUT}.tbi" ]; then
  echo "TABIXING INPUT DATA"
  tabix ${INPUT}
fi

#REFERENCE PATH
if [ ! -e "${REFERENCE}.vcf.gz" ]; then
  echo "CREATION OF REFERENCE IN VCF FORMAT"
  bcftools view ${REFERENCE}.bcf.gz -Oz -o ${REFERENCE}.vcf.gz && bcftools index ${REFERENCE}.vcf.gz
fi
if [ ! -e "${REFERENCE}.bcf.gz" ]; then
  echo "CREATION OF REFERENCE IN BCF FORMAT" 
  bcftools view ${REFERENCE}.vcf.gz -Ob -o ${REFERENCE}.bcf.gz && bcftools index ${REFERENCE}.bcf.gz
fi


ref_eagle2=${REFERENCE}.bcf.gz
ref_shapeit4=${REFERENCE}.bcf.gz


if [ -f "${REFERENCE}.bref3" ]; then
  ref_beagle5=${REFERENCE}.bref3
  
else
  echo "missing bref3 format for beagle, creating...."
  set -e
  $bref3 ${REFERENCE}.vcf.gz > ${REFERENCE}.bref3
  echo "Created ${REFERENCE}.bref3"
  ref_beagle5=${REFERENCE}.bref3
fi

#CHR
CHR=${CHROMOSOME}
#OUTPUT PATH
if [ ! -d ${base_output} ]; then
  mkdir -p ${base_output};
  if [ ! -d ${base_output}/phasing ]; then
    mkdir -p ${base_output}/phasing;
    if [ ! -d ${base_output}/phasing/eagle ]; then
      mkdir -p ${base_output}/phasing/eagle;
    fi
    if [ ! -d ${base_output}/phasing/shapeit ]; then
      mkdir -p ${base_output}/phasing/shapeit;
    fi
    if [ ! -d ${base_output}/phasing/beagle_phase ]; then
      mkdir -p ${base_output}/phasing/beagle_phase;
    fi
  fi
else
  echo "[results] folder already exists! "
  set -e
fi


output_eagle2=${base_output}/phasing/eagle/${OUTPUT}
output_shapeit4=${base_output}/phasing/shapeit/${OUTPUT}
output_beagle5_phase=${base_output}/phasing/beagle_phase/${OUTPUT}


###################### PHASING ########################
#######################################################

if [ $EAGLE == "true" ]; then
  echo "EAGLE with REF"
  $time "EAGLE with REF MEMORY= %M" $eagle2 \
      --vcfRef=${ref_eagle2} \
      --vcfTarget=${INPUT} \
      --geneticMapFile=${map_eagle2} \
      --chrom=${CHR} \
      --outPrefix=${output_eagle2}_chr${CHR}_eaglePhasing_REF \
      --numThreads=${THREADS} \
      2>&1 | tee ${output_eagle2}_chr${CHR}_eaglePhasing_REF.log
  
  echo "EAGLE NO REF"
  $time "EAGLE NO REF MEMORY= %M" $eagle2 \
      --vcf=${INPUT} \
      --geneticMapFile=${map_eagle2} \
      --chrom=${CHR} \
      --outPrefix=${output_eagle2}_chr${CHR}_eaglePhasing_noREF \
      --numThreads=${THREADS} \
      2>&1 | tee ${output_eagle2}_chr${CHR}_eaglePhasing_noREF.log
  
fi
if [ $SHAPEIT == "true" ]; then
  echo "SHAPEIT NO REF"
  $time "SHAPEIT NO REF MEMORY= %M" $shapeit4 \
      --input ${INPUT} \
      --map ${map_shapeit4} \
      --region ${CHR} \
      --output ${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF.vcf.gz \
      --thread ${THREADS} \
      --log ${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF.log
  
  echo "SHAPEIT with REF"
  $time "SHAPEIT with REF MEMORY= %M" $shapeit4 \
      --input ${INPUT} \
      --map ${map_shapeit4} \
      --region ${CHR} \
      --reference ${ref_shapeit4} \
      --output ${output_shapeit4}_chr${CHR}_shapeitPhasing_REF.vcf.gz \
      --thread ${THREADS} \
      --log ${output_shapeit4}_chr${CHR}_shapeitPhasing_REF.log
  
fi
if [ $BEAGLE_PHASE == "true" ]; then
  echo "BEAGLE with REF"
  $time "BEAGLE with REF MEMORY= %M" $beagle5 \
      gt=${INPUT} \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_phase}_chr${CHR}_beaglePhasing_REF \
      impute=false 
  
  echo "BEAGLE NO REF"
  $time "BEAGLE NO REF MEMORY= %M" $beagle5 \
      gt=${INPUT} \
      map=${map_beagle5} \
      out=${output_beagle5_phase}_chr${CHR}_beaglePhasing_noREF \
      impute=false 
  
fi

for i in `ls ${base_output}/phasing/*/${OUTPUT}*.vcf.gz`; do
  tabix ${i} -f
done

echo "Accuracy ~°~°~° Accuracy ~°~°~° Accuracy ~°~°~° Accuracy ~°~°~°"
######## Accuracy ######## Accuracy ########### Accuracy ########### Accuracy ##############
for i in `ls ${base_output}/phasing/*/${OUTPUT}*.vcf.gz`; do 
    zcat $i | grep -v "##" | cut -f 3,10- | sed 's/0|0/0/g'|sed 's/1|0/1/g'|sed 's/0|1/2/g'|sed 's/1|1/3/g' > ${i}.numpy
done

