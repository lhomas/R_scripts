#!/bin/bash


# USAGE: /path/to/bamStatGen.sh /path/to/file.bam sample chromosome start_interval end_interval
# This file takes a bam or cram file as input and output coverage and mapping quality statistics over a specified region of the genome. 

samtools view ${1} ${3}:${4}-${5} > ${2}.interval.bam



mapq_end=$((${5}-10000)) 

echo ${mapq_end}

seq -f %1.0f ${4} 10000 ${mapq_end} | \
while read r
do
	samtools view ${1} ${3}:${r}-$((${r}+10000)) >> ${2}_mapq.txt
	samtools view ${1} ${3}:${r}-$((${r}+10000)) | \
	awk '{sum+=$5} END {print NR ? sum/NR : 0}'
done > ${2}.MAPQ.${3}_${4}_${5}_10000.txt


echo "COVERAGE"
cov_end=$((${5}-1000))

seq -f %1.0f ${4} 1000 ${cov_end} | \
while read r
do 
	samtools depth -r ${3}:${r}-$((${r}+1000)) ${1} >> ${2}_cov.txt
	samtools depth -r ${3}:${r}-$((${r}+1000)) ${1} | \
	awk '{sum+=$3} END {print NR ? sum/NR : 0}'
done > ${2}.cov.${3}_${4}_${5}_100.txt

