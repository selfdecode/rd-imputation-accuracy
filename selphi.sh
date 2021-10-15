#!/bin/bash
if [ -z "$*" ]; then 
    echo "
author: SELFDECODE
contact email : adriano@selfdecode.com

SelfDecode software to analyze multiple phasing and imputation softwares

Parameters required:
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
              chomosome to analyze allowed [1-22 X]. NO chr prefix."
      echo"
Optional Parameters:      
"
      echo "      -ibeagle|--imp_beagle
              skip Beagle imputation"
      echo "      -pbeagle|--phase_beagle
              skip Beagle haplotype estimation"
      echo "      -impute5|--impute5
              skip Impute imputation"
      echo "      -shapeit|--shapeit
              skip ShapeIT Phasing"
      echo "      --minimac|--minimac
              skip Minimac imputation"
      echo "      -eagle|--eagle
              skip Eagle phasing"
      echo "      -bigref|--BIGREF
              use this option if you get memory allocate error during accuracy evaluation"

echo "
[base] Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>
[skip] Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -ibeagle x -impute5 x
[memo] Usage: ./PIscore.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -bigref x

"
    exit 1
fi

BEAGLE_IMP=true
BEAGLE_PHASE=true
IMPUTE5=true
MINIMAC=true
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
    -ibeagle|--imp_beagle)
      BEAGLE_IMP="false"
      shift 
      shift
      ;;
    -pbeagle|--phase_beagle)
      BEAGLE_PHASE="false"
      shift 
      shift
      ;;
    -impute5|--impute5)
      IMPUTE5="false"
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
    -minimac|--minimac)
      MINIMAC="false"
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
              chomosome to analyze allowed [1-22 X]. NO chr prefix."
      echo "      -ibeagle|--imp_beagle
              skip Beagle imputation"
      echo "      -pbeagle|--phase_beagle
              skip Beagle haplotype estimation"
      echo "      -impute5|--impute5
              skip Impute imputation"
      echo "      -shapeit|--shapeit
              skip ShapeIT Phasing"
      echo "      --minimac|--minimac
              skip Minimac imputation"
      echo "      -eagle|--eagle
              skip Eagle phasing"
      echo "      -bigref|--BIGREF
              use this option if you get memory allocate error during accuracy evaluation"

echo "
[base] Usage: ./ImputePermute.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20>
[skip] Usage: ./ImputePermute.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -ibeagle x -impute5 x
[memo] Usage: ./ImputePermute.sh -i <input.vcf.gz> -r <ref_file> -t <4> -o <output_name> -c <20> -bigref x

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
echo "BEAGLE IMP    : ${BEAGLE_IMP}"
echo "BEAGLE PHASE  : ${BEAGLE_PHASE}"
echo "EAGLE         : ${EAGLE}"
echo "IMPUTE5       : ${IMPUTE5}"
echo "SHAPEIT       : ${SHAPEIT}"
echo "MINIMAC       : ${MINIMAC}"
echo "BIGREF        : ${BIGREF}"

#######################################################
#######################################################
#SOFTWARE PATH (EDITABLE PART)
time="/usr/bin/time -f"
eagle2="/home/ec2-user/adriano/imputation/phase2/software/eagle2.4.1/Eagle_v2.4.1/eagle"
shapeit4="/home/ec2-user/adriano/imputation/phase2/software/shapeit4/shapeit4-4.2.1/bin/shapeit4.2"
beagle5="java -Xmx8g -jar /home/ec2-user/adriano/imputation/phase2/software/beagle5.2/beagle.29May21.d6d.jar"
bref3="java -Xmx8g -jar /home/ec2-user/adriano/imputation/phase2/software/beagle5.2/bref3.29May21.d6d.jar"
imp5Converter="/home/ec2-user/adriano/imputation/phase2/software/impute5/impute5_v1.1.5/imp5Converter_1.1.5_static"
miniConverter="/home/ec2-user/adriano/imputation/phase2/software/Minimac3Executable/bin/Minimac3"
impute5="/home/ec2-user/adriano/imputation/phase2/software/impute5/impute5_v1.1.5/impute5_1.1.5_static"
minimac4="/home/ec2-user/adriano/imputation/phase2/software/minimac4/Minimac4/build/minimac4"
simpy="/home/ec2-user/adriano/git/rd-imputation-accuracy/Simpy.py"
imputation_accuracy="/home/ec2-user/adriano/git/rd-imputation-accuracy/imputation_accuracy.sh"
#GENETIC RECOMBINATIO MAP PATH
map_eagle2="/home/ec2-user/adriano/imputation/phase2/software/eagle2.4.1/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
map_beagle5="/home/ec2-user/adriano/imputation/phase2/genetic_map/plink.chr${CHROMOSOME}.chr.GRCh38.map"
map_shapeit4="/home/ec2-user/adriano/imputation/phase2/genetic_map/chr${CHROMOSOME}.b38.gmap.gz"
map_impute5=$map_shapeit4
#WGS PATH for accuracy
wgs_subset_chr="/home/ec2-user/adriano/imputation/phase2/reference_panel/ref_30x/chr20.reference_panel.30x.hg38.190samples.vcf.gz"
wgs_subset="/home/ec2-user/adriano/imputation/phase2/reference_panel/ref_30x/20.reference_panel.30x.hg38.190samples.vcf.gz"
bwgs_subset_chr="/home/ec2-user/adriano/imputation/phase2/reference_panel/ref_30x/chr20.reference_panel.30x.hg38.190samples.bcf.gz"
bwgs_subset="/home/ec2-user/adriano/imputation/phase2/reference_panel/ref_30x/20.reference_panel.30x.hg38.190samples.bcf.gz"
#######################################################
#######################################################

if [ ! -e "${wgs_subset}" ]; then
  echo "CREATION OF WGS SUBSET IN VCF FORMAT"
  echo "chr${CHROMOSOME} ${CHROMOSOME}" > chr${CHROMOSOME}_${CHROMOSOME}.txt
  bcftools annotate --rename-chrs chr${CHROMOSOME}_${CHROMOSOME}.txt ${wgs_subset_chr} \
    -Oz -o ${wgs_subset} \
    --thread ${THREADS} && tabix ${output_name} -f
  bcftools view ${wgs_subset} -Ob -o ${bwgs_subset} && bcftools index ${bwgs_subset}
  rm chr${CHROMOSOME}_${CHROMOSOME}.txt
fi

#REFERENCE PATH
if [ ! -e "${REFERENCE}.vcf.gz" ]; then
  bcftools view ${REFERENCE}.bcf.gz -Ov -o ${REFERENCE}.vcf.gz && bcftools index ${REFERENCE}.vcf.gz
fi
if [ ! -e "${REFERENCE}.bcf.gz" ]; then
  bcftools view ${REFERENCE}.vcf.gz -Ov -o ${REFERENCE}.bcf.gz && bcftools index ${REFERENCE}.bcf.gz
fi


ref_eagle2=${REFERENCE}.bcf.gz
ref_shapeit4=${REFERENCE}.bcf.gz


if [ -f "${REFERENCE}.bref3" ]; then
  ref_beagle5=${REFERENCE}.bref3
else
  echo "missing bref3 format for beagle, creating...."
  $bref3 ${REFERENCE}.vcf.gz > ${REFERENCE}.bref3
  echo "Created ${REFERENCE}.bref3"
  ref_beagle5=${REFERENCE}.bref3
fi

if [ -f "${REFERENCE}.imp5" ]; then
  ref_impute5=${REFERENCE}.imp5
else
  echo "missing imp5 format for impute5, creating...."
  $imp5Converter \
    --h ${REFERENCE}.vcf.gz \
    --r chr${CHROMOSOME} \
    --o ${REFERENCE}.imp5
  echo "Created ${REFERENCE}.imp5"
  ref_impute5=${REFERENCE}.imp5
fi

if [ -f "${REFERENCE}.m3vcf.gz" ]; then
  ref_minimac4=${REFERENCE}.m3vcf.gz
else
  echo "missing m3vcf format for minimac4, creating....
  Minimac accept only vcf format with conting name [1-22] X"
  if [ -f "${REFERENCE}.vcf.gz" ]; then
    echo "vcf format exist...checking the CONTIG format"
    STR=`bcftools query -f '%CHROM\n' ${REFERENCE}.vcf.gz | head -1`
    SUB='chr'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo "CONTIG name with chr"
      echo "chr${CHROMOSOME} ${CHROMOSOME}" > chr${CHROMOSOME}_${CHROMOSOME}.txt
      bcftools annotate --rename-chrs chr${CHROMOSOME}_${CHROMOSOME}.txt ${REFERENCE}.vcf.gz -Oz -o ${REFERENCE}_minimac_tmp.vcf.gz --threads ${THREADS}
      rm chr${CHROMOSOME}_${CHROMOSOME}.txt
    fi 
  else
    echo "converting BCF in VCF for minimac"
    bcftools view ${REFERENCE}.bcf.gz -Oz -o ${REFERENCE}.vcf.gz
    STR=`bcftools query -f '%CHROM\n' ${REFERENCE}.vcf.gz | head -1`
    SUB='chr'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo "CONTIG name with chr...converting to [1-22] X"
      echo "chr${CHROMOSOME} ${CHROMOSOME}" > chr${CHROMOSOME}_${CHROMOSOME}.txt
      bcftools annotate --rename-chrs chr${CHROMOSOME}_${CHROMOSOME}.txt ${REFERENCE}.vcf.gz -Oz -o ${REFERENCE}_minimac_tmp.vcf.gz --threads ${THREADS}
      rm chr${CHROMOSOME}_${CHROMOSOME}.txt
    fi 
  fi
  $miniConverter \
    --refHaps ${REFERENCE}_minimac_tmp.vcf.gz \
    --processReference \
    --prefix ${REFERENCE}
  echo "Created ${REFERENCE}.m3vcf.gz"
  ref_minimac4=${REFERENCE}.m3vcf.gz
fi
#CHR
CHR=${CHROMOSOME}
#OUTPUT PATH
base_output=/home/ec2-user/adriano/imputation/phase2/results
output_eagle2=${base_output}/phasing/eagle/${OUTPUT}
output_shapeit4=${base_output}/phasing/shapeit/${OUTPUT}
output_beagle5_phase=${base_output}/phasing/beagle_phase/${OUTPUT}
output_beagle5_imp=${base_output}/imputation/beagle_imp/${OUTPUT}
output_impute5=${base_output}/imputation/impute5/${OUTPUT}
output_minimac4=${base_output}/imputation/minimac/${OUTPUT}

#######################################################
#######################################################


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
      --region chr${CHR} \
      --output ${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF.vcf.gz \
      --thread ${THREADS} \
      --log ${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF.log
  
  echo "SHAPEIT with REF"
  $time "SHAPEIT with REF MEMORY= %M" $shapeit4 \
      --input ${INPUT} \
      --map ${map_shapeit4} \
      --region chr${CHR} \
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


#################### IMPUTATION #######################
#######################################################

#################### BEAGLE #######################
if [ $BEAGLE_IMP == "true" ]; then
  echo "BEAGLE-noREF-BEAGLE"
  $time "BEAGLE-noREF-BEAGLE MEMORY= %M" $beagle5 \
      gt=${output_beagle5_phase}_chr${CHR}_beaglePhasing_noREF.vcf.gz \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_imp}_chr${CHR}_beaglePhasing_noREF_beagleImputed \
      gp=true \
      ap=true
  
  echo "BEAGLE-REF-BEAGLE"
  $time "BEAGLE-REF-BEAGLE MEMORY= %M" $beagle5 \
      gt=${output_beagle5_phase}_chr${CHR}_beaglePhasing_REF.vcf.gz \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_imp}_chr${CHR}_beaglePhasing_REF_beagleImputed \
      gp=true \
      ap=true
  
  echo "EAGLE-noREF-BEAGLE"
  $time "EAGLE-noREF-BEAGLE MEMORY= %M" $beagle5 \
      gt=${output_eagle2}_chr${CHR}_eaglePhasing_noREF.vcf.gz \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_imp}_chr${CHR}_eaglePhasing_noREF_beagleImputed \
      gp=true \
      ap=true
  
  echo "EAGLE-REF-BEAGLE"
  $time "EAGLE-REF-BEAGLE MEMORY= %M" $beagle5 \
      gt=${output_eagle2}_chr${CHR}_eaglePhasing_REF.vcf.gz \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_imp}_chr${CHR}_eaglePhasing_REF_beagleImputed \
      gp=true \
      ap=true
  
  echo "SHAPEIT-noREF-BEAGLE"
  $time "SHAPEIT-noREF-BEAGLE MEMORY= %M" $beagle5 \
      gt=${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF.vcf.gz \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_imp}_chr${CHR}_shapeitPhasing_noREF_beagleImputed \
      gp=true \
      ap=true
  
  echo "SHAPEIT-REF-BEAGLE"
  $time "SHAPEIT-REF-BEAGLE MEMORY= %M" $beagle5 \
      gt=${output_shapeit4}_chr${CHR}_shapeitPhasing_REF.vcf.gz \
      map=${map_beagle5} \
      ref=${ref_beagle5} \
      out=${output_beagle5_imp}_chr${CHR}_shapeitPhasing_REF_beagleImputed \
      gp=true \
      ap=true
fi

if [ $IMPUTE5 == "true" ]; then
  #################### IMPUTE5 #######################
  echo "EAGLE-REF-IMPUTE5"
  $time "EAGLE-REF-IMPUTE5 MEMORY= %M" $impute5 \
    --h ${ref_impute5} \
    --m ${map_impute5} \
    --g ${output_eagle2}_chr${CHR}_eaglePhasing_REF.vcf.gz \
    --r chr${CHR} \
    --o ${output_impute5}_chr${CHR}_eaglePhasing_REF_impute5Imputed.vcf.gz \
    --threads ${THREADS} \
    --l ${output_impute5}_chr${CHR}_eaglePhasing_REF_impute5Imputed.log
  
  echo "EAGLE-noREF-IMPUTE5"
  $time "EAGLE-noREF-IMPUTE5 MEMORY= %M" $impute5 \
    --h ${ref_impute5} \
    --m ${map_impute5} \
    --g ${output_eagle2}_chr${CHR}_eaglePhasing_noREF.vcf.gz \
    --r chr${CHR} \
    --o ${output_impute5}_chr${CHR}_eaglePhasing_noREF_impute5Imputed.vcf.gz \
    --threads ${THREADS} \
    --l ${output_impute5}_chr${CHR}_eaglePhasing_noREF_impute5Imputed.log
  
  echo "SHAPEIT-REF-IMPUTE5"
  $time "SHAPEIT-REF-IMPUTE5 MEMORY= %M" $impute5 \
    --h ${ref_impute5} \
    --m ${map_impute5} \
    --g ${output_shapeit4}_chr${CHR}_shapeitPhasing_REF.vcf.gz \
    --r chr${CHR} \
    --o ${output_impute5}_chr${CHR}_shapeitPhasing_REF_impute5Imputed.vcf.gz \
    --threads ${THREADS} \
    --l ${output_impute5}_chr${CHR}_shapeitPhasing_REF_impute5Imputed.log
  
  echo "SHAPEIT-noREF-IMPUTE5"
  $time "SHAPEIT-noREF-IMPUTE5 MEMORY= %M" $impute5 \
    --h ${ref_impute5} \
    --m ${map_impute5} \
    --g ${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF.vcf.gz \
    --r chr${CHR} \
    --o ${output_impute5}_chr${CHR}_shapeitPhasing_noREF_impute5Imputed.vcf.gz \
    --threads ${THREADS} \
    --l ${output_impute5}_chr${CHR}_shapeitPhasing_noREF_impute5Imputed.log
  
  echo "BEAGLE-REF-IMPUTE5"
  $time "BEAGLE-REF-IMPUTE5 MEMORY= %M" $impute5 \
    --h ${ref_impute5} \
    --m ${map_impute5} \
    --g ${output_beagle5_phase}_chr${CHR}_beaglePhasing_REF.vcf.gz \
    --r chr${CHR} \
    --o ${output_impute5}_chr${CHR}_beaglePhasing_REF_impute5Imputed.vcf.gz \
    --threads ${THREADS} \
    --l ${output_impute5}_chr${CHR}_beaglePhasing_REF_impute5Imputed.log
  
  echo "BEAGLE-noREF-IMPUTE5"
  $time "BEAGLE-noREF-IMPUTE5 MEMORY= %M" $impute5 \
    --h ${ref_impute5} \
    --m ${map_impute5} \
    --g ${output_beagle5_phase}_chr${CHR}_beaglePhasing_noREF.vcf.gz \
    --r chr${CHR} \
    --o ${output_impute5}_chr${CHR}_beaglePhasing_noREF_impute5Imputed.vcf.gz \
    --threads ${THREADS} \
    --l ${output_impute5}_chr${CHR}_beaglePhasing_noREF_impute5Imputed.log
fi

#################### MINIMAC #######################
#preparing data for minimac that only handle data without chr prefix
if [ $MINIMAC == "true" ]; then
  echo "chr${CHROMOSOME} ${CHROMOSOME}" > chr${CHROMOSOME}_${CHROMOSOME}.txt
  for i in `ls ${base_output}/phasing/*/${OUTPUT}*.vcf.gz`; do
    output_name=${i//.vcf.gz/_minimac.vcf.gz}
    bcftools annotate --rename-chrs chr${CHROMOSOME}_${CHROMOSOME}.txt $i \
      -Oz -o ${output_name} \
      --thread ${THREADS} && tabix ${output_name} -f
  done
  rm chr${CHROMOSOME}_${CHROMOSOME}.txt
  echo "EAGLE-REF-MINIMAC"
  $time "EAGLE-REF-MINIMAC MEMORY= %M" $minimac4 \
    --refHaps ${ref_minimac4} \
    --haps ${output_eagle2}_chr${CHR}_eaglePhasing_REF_minimac.vcf.gz \
    --prefix ${output_minimac4}_chr${CHR}_eaglePhasing_REF_minimacImputed \
    --cpus ${THREADS} \
    --log 
  
  echo "EAGLE-noREF-MINIMAC"
  $time "EAGLE-noREF-MINIMAC MEMORY= %M" $minimac4 \
    --refHaps ${ref_minimac4} \
    --haps ${output_eagle2}_chr${CHR}_eaglePhasing_noREF_minimac.vcf.gz \
    --prefix ${output_minimac4}_chr${CHR}_eaglePhasing_noREF_minimacImputed \
    --cpus ${THREADS} \
    --log 
  
  echo "SHAPEIT-REF-MINIMAC"
  $time "SHAPEIT-REF-MINIMAC MEMORY= %M" $minimac4 \
    --refHaps ${ref_minimac4} \
    --haps ${output_shapeit4}_chr${CHR}_shapeitPhasing_REF_minimac.vcf.gz \
    --prefix ${output_minimac4}_chr${CHR}_shapeitPhasing_REF_minimacImputed \
    --cpus ${THREADS} \
    --log 
  
  echo "SHAPEIT-noREF-MINIMAC"
  $time "SHAPEIT-noREF-MINIMAC MEMORY= %M" $minimac4 \
    --refHaps ${ref_minimac4} \
    --haps ${output_shapeit4}_chr${CHR}_shapeitPhasing_noREF_minimac.vcf.gz \
    --prefix ${output_minimac4}_chr${CHR}_shapeitPhasing_noREF_minimacImputed \
    --cpus ${THREADS} \
    --log 
  
  echo "BEAGLE-REF-MINIMAC"
  $time "BEAGLE-REF-MINIMAC MEMORY= %M" $minimac4 \
    --refHaps ${ref_minimac4} \
    --haps ${output_beagle5_phase}_chr${CHR}_beaglePhasing_REF_minimac.vcf.gz \
    --prefix ${output_minimac4}_chr${CHR}_beaglePhasing_REF_minimacImputed \
    --cpus ${THREADS} \
    --log 
  
  echo "BEAGLE-noREF-MINIMAC"
  $time "BEAGLE-noREF-MINIMAC MEMORY= %M" $minimac4 \
    --refHaps ${ref_minimac4} \
    --haps ${output_beagle5_phase}_chr${CHR}_beaglePhasing_noREF_minimac.vcf.gz \
    --prefix ${output_minimac4}_chr${CHR}_beaglePhasing_noREF_minimacImputed \
    --cpus ${THREADS} \
    --log 
fi

echo "Accuracy ~°~°~° Accuracy ~°~°~° Accuracy ~°~°~° Accuracy ~°~°~°"
######## Accuracy ######## Accuracy ########### Accuracy ########### Accuracy ##############
if [ $BIGREF == "true" ]; then
  for i in `ls ${base_output}/imputation/*/${OUTPUT}*.vcf.gz`; do   
    echo $i
    if [ ! -e "${i}.tbi" ]; then
      tabix $i -f
    fi
    STR=$i
    SUB='minimac'
    if [[ "$STR" != *"$SUB"* ]]; then
      $imputation_accuracy \
        -i $i \
        -w $wgs_subset_chr \
        -bw $bwgs_subset_chr \
        -t ${THREADS}
    else
      $imputation_accuracy \
        -i $i \
        -w $wgs_subset \
        -bw $bwgs_subset \
        -t ${THREADS}
    fi
  done
else
  for i in `ls ${base_output}/imputation/*/${OUTPUT}*.vcf.gz`; do   
    echo $i
    if [ ! -e "${i}.tbi" ]; then
      tabix $i -f
    fi
    STR=$i
    SUB='minimac'
    if [[ "$STR" != *"$SUB"* ]]; then
      python3 $simpy \
        --imputed $i \
        --wgs $wgs_subset_chr \
        --bwgs $bwgs_subset_chr
    else
      python3 $simpy \
        --imputed $i \
        --wgs $wgs_subset \
        --bwgs $bwgs_subset
    fi
  done
fi
