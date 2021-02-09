#!/bin/bash
#set -o errexit
#set -o pipefail

# Help message --------------------------------------------------------------------------
# Require version 4 of bash which have associative arrays
if [[ $BASH_VERSINFO < 4 ]]; 
    then echo "bash version is out of date, requires at least version 4"; 
fi

if [ $# -eq 0 ]; then
echo >&2 "
$(basename $0) - Merge many bam files and then generate bigwigs.
USAGE: bash $(basename $0) -p <cores> -c <merge configuration> -o <output directory>

-p  Number of cores to use per merge process            [default = 4]
-c  Configuration file for the merge                    [required]
-o  Output directory for merged bam and bigwigs         [required]
"
exit 1
fi

# Parse arguments -----------------------------------------------------------------------
p=4
while getopts "p:c:o:" op
do
    case "$op" in
        p)  p="$OPTARG";;
        c)  c="$OPTARG";;
        o)  o="$OPTARG";;
        \?) exit 1;;
    esac
done

if [[ ! -e $o ]]; then
    mkdir -p $o
fi

log="mergebams.log"
if [[ -e $log ]]; then
    rm $log
fi
touch $log

# Associate target file with starting files ---------------------------------------------
echo -e "Associating starting files with target files\n"
declare -A merge_dict
while IFS=$'\t' read starting_file target_file
do
    # skip the header
    if [[ $starting_file == "bam" ]]; then
        continue
    fi
    # echo $target_file, $starting_file
    merge_dict+=( [$target_file]=" $starting_file" )
done < $c

# merge_dict key is target file, and values are starting files separated by a space.
echo -e "target_file \t starting_file"
for key in ${!merge_dict[@]}; do
    echo -e "${key}:${merge_dict[${key}]}"
done

# merge and index files -----------------------------------------------------------------
echo -e "\nMerging and sorting bam files\n"
for target in ${!merge_dict[@]}; do
    in_files=${merge_dict[${target}]}
    out_file=$o/$target
    if test -f "$out_file"; then
        echo "$out_file exists."
        continue
    fi
    echo "samtools merge -@ $p - $in_files | samtools sort -@ $p - -o $out_file &" >> $log
    samtools merge -@ $p - $in_files | samtools sort -@ $p - -o $out_file &
done
wait

echo -e "Indexing bam files\n"
for target in ${!merge_dict[@]}; do
    in_file=$o/$target
    out_file=$o/$target.bai
    if test -f "$out_file"; then
        echo "$out_file exists."
        continue
    fi
    echo "samtools index -@ $p $in_file &" >> $log
    samtools index -@ $p $in_file &
done
wait

# make bigwigs --------------------------------------------------------------------------
echo -e "Making bigwigs\n"
for target in ${!merge_dict[@]}; do
    in_file=$o/$target
    out_file=$o/$(basename $target).bw
    if test -f "$out_file"; then
        echo "$out_file" exists.
        continue
    fi
    echo "bamCoverage -b $in_file -o $out_file -p $p --binSize 10 --smoothLength 50 --normalizeUsing CPM &" >> $log
    bamCoverage -b $in_file -o $out_file -p $p --binSize 10 --smoothLength 50 --normalizeUsing CPM &
done

echo -e "Waiting for bigwigs to generate\n"
wait
echo -e "Bam files merged and bigwigs made. Program is finished.\n"