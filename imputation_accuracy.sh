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

# Usage of the progress bar:
# Source this script
# enable_trapping <- optional to clean up properly if user presses ctrl-c
# setup_scroll_area <- create empty progress bar
# draw_progress_bar 10 <- advance progress bar
# draw_progress_bar 40 <- advance progress bar
# block_progress_bar 45 <- turns the progress bar yellow to indicate some action is requested from the user
# draw_progress_bar 90 <- advance progress bar
# destroy_scroll_area <- remove progress bar

# Constants
CODE_SAVE_CURSOR="\033[s"
CODE_RESTORE_CURSOR="\033[u"
CODE_CURSOR_IN_SCROLL_AREA="\033[1A"
COLOR_FG="\e[30m"
COLOR_BG="\e[42m"
COLOR_BG_BLOCKED="\e[43m"
RESTORE_FG="\e[39m"
RESTORE_BG="\e[49m"

# Variables
PROGRESS_BLOCKED="false"
TRAPPING_ENABLED="false"
TRAP_SET="false"

CURRENT_NR_LINES=0

setup_scroll_area() {
    # If trapping is enabled, we will want to activate it whenever we setup the scroll area and remove it when we break the scroll area
    if [ "$TRAPPING_ENABLED" = "true" ]; then
        trap_on_interrupt
    fi

    lines=$(tput lines)
    CURRENT_NR_LINES=$lines
    let lines=$lines-1
    # Scroll down a bit to avoid visual glitch when the screen area shrinks by one row
    echo -en "\n"

    # Save cursor
    echo -en "$CODE_SAVE_CURSOR"
    # Set scroll region (this will place the cursor in the top left)
    echo -en "\033[0;${lines}r"

    # Restore cursor but ensure its inside the scrolling area
    echo -en "$CODE_RESTORE_CURSOR"
    echo -en "$CODE_CURSOR_IN_SCROLL_AREA"

    # Start empty progress bar
    draw_progress_bar 0
}

destroy_scroll_area() {
    lines=$(tput lines)
    # Save cursor
    echo -en "$CODE_SAVE_CURSOR"
    # Set scroll region (this will place the cursor in the top left)
    echo -en "\033[0;${lines}r"

    # Restore cursor but ensure its inside the scrolling area
    echo -en "$CODE_RESTORE_CURSOR"
    echo -en "$CODE_CURSOR_IN_SCROLL_AREA"

    # We are done so clear the scroll bar
    clear_progress_bar

    # Scroll down a bit to avoid visual glitch when the screen area grows by one row
    echo -en "\n\n"

    # Once the scroll area is cleared, we want to remove any trap previously set. Otherwise, ctrl+c will exit our shell
    if [ "$TRAP_SET" = "true" ]; then
        trap - INT
    fi
}

draw_progress_bar() {
    percentage=$1
    lines=$(tput lines)
    let lines=$lines

    # Check if the window has been resized. If so, reset the scroll area
    if [ "$lines" -ne "$CURRENT_NR_LINES" ]; then
        setup_scroll_area
    fi

    # Save cursor
    echo -en "$CODE_SAVE_CURSOR"

    # Move cursor position to last row
    echo -en "\033[${lines};0f"

    # Clear progress bar
    tput el

    # Draw progress bar
    PROGRESS_BLOCKED="false"
    print_bar_text $percentage

    # Restore cursor position
    echo -en "$CODE_RESTORE_CURSOR"
}

block_progress_bar() {
    percentage=$1
    lines=$(tput lines)
    let lines=$lines
    # Save cursor
    echo -en "$CODE_SAVE_CURSOR"

    # Move cursor position to last row
    echo -en "\033[${lines};0f"

    # Clear progress bar
    tput el

    # Draw progress bar
    PROGRESS_BLOCKED="true"
    print_bar_text $percentage

    # Restore cursor position
    echo -en "$CODE_RESTORE_CURSOR"
}

clear_progress_bar() {
    lines=$(tput lines)
    let lines=$lines
    # Save cursor
    echo -en "$CODE_SAVE_CURSOR"

    # Move cursor position to last row
    echo -en "\033[${lines};0f"

    # clear progress bar
    tput el

    # Restore cursor position
    echo -en "$CODE_RESTORE_CURSOR"
}

print_bar_text() {
    local percentage=$1
    local cols=$(tput cols)
    let bar_size=$cols-19

    local color="${COLOR_FG}${COLOR_BG}"
    if [ "$PROGRESS_BLOCKED" = "true" ]; then
        color="${COLOR_FG}${COLOR_BG_BLOCKED}"
    fi

    # Prepare progress bar
    let complete_size=($bar_size*$percentage)/$n_file
    let remainder_size=$bar_size-$complete_size
    progress_bar=$(echo -ne "["; echo -en "${color}"; printf_new "#" $complete_size; echo -en "${RESTORE_FG}${RESTORE_BG}"; printf_new "." $remainder_size; echo -ne "]");

    # Print progress bar
    echo -ne " Progress ${percentage}/$n_file ${progress_bar}"
}

enable_trapping() {
    TRAPPING_ENABLED="true"
}

trap_on_interrupt() {
    # If this function is called, we setup an interrupt handler to cleanup the progress bar
    TRAP_SET="true"
    trap cleanup_on_interrupt INT
}

cleanup_on_interrupt() {
    destroy_scroll_area
    exit
}

printf_new() {
    str=$1
    num=$2
    v=$(printf "%-${num}s" "$str")
    echo -ne "${v// /$str}"
}


if [ $skip_chunks == "False" ]; then 
  #grab the header
  zcat ${IMPUTED} | head -10000 | grep '^#' > header
  #grab the non header lines
  zgrep -v '^#' ${IMPUTED} > variants
  #split into chunks with 1000 lines
  split -l 100000 variants
  #reattach the header to each and clean up
  for i in x*;do 
    cat header $i | bgzip > $i.vcf.gz
    tabix -p vcf $i.vcf.gz
    rm -f $i;
  done
  echo "Splitted Imputed file in chuncks of [100k]"
  rm -f header variants
  for i in x*.vcf.gz; do 
    filename=$(basename -- "$i")
    out_name="${filename%%.*}"
    bcftools view $i -Ob -o $out_name.bcf.gz --threads 4
    bcftools index $out_name.bcf.gz;
  done
  echo "BCF Imputed files Created"
fi

if [ $drop_analysis == "False" ]; then 
  n_file=`ls x*.vcf.gz | wc -l`
  counter=1; 
  #setup_scroll_area
  for splitted in x*.vcf.gz; do
      #printf "processing $((counter++)) of $n_file files\n"
      draw_progress_bar $((counter++))
      filename=$(basename -- "$splitted")
      file_name="${filename%%.*}"
      max_bg_procs $THREADS
      python3 /home/ec2-user/adriano/git/imputation_accuracy_calculator/Simpy.py --wgs ${WGS} --imputed ${file_name}.vcf.gz --bwgs ${BWGS} --bimputed ${file_name}.bcf.gz &
  done
  wait
  destroy_scroll_area
  echo "Joining files..."
  filename=$(basename -- "$IMPUTED")
  out_name="${filename%%.*}"
  awk 'FNR>1 || NR==1' *_per_sample_results* | bgzip > ${out_name}_per_sample_results.tsv.gz
  awk 'FNR>1 || NR==1' *_per_variant_results* | bgzip > ${out_name}_per_variant_results.tsv.gz
  #delete tmp files
  echo "Deleting tmp files..."
  rm x*.vcf.gz*
  rm x*.bcf.gz*
  rm x*.vcf_per_variant_results*
  rm x*.vcf_per_sample_results*
fi

python3 /home/ec2-user/adriano/git/imputation_accuracy_calculator/rebuild_metrics.py --s ${out_name}_per_sample_results.tsv.gz

awk 'FNR>1 || NR==1' *_ImputationAccuracy.txt  > ${out_name}_per_sample_results.tsv
rm *_ImputationAccuracy.txt
rm ${out_name}_per_sample_results.tsv.gz
echo "Output-1    : ${out_name}_per_sample_results.tsv"
echo "Output-2    : ${out_name}_per_variant_results.tsv.gz"