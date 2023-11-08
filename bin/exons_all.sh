#!/bin/bash

dir=`dirname $0`
bam=gene_design_uniq.bam
bam=gene_qual_design_uniq.bam
samtools view -H $bam

if ! [ -e gene_exons.bed ]; then
  echo "$dir/exons.pl $dir/../ref/hg38_exons.bed > gene_exons.bed"
  exit -1
fi

for g in `awk '{print $4}' gene_exons.bed | uniq`; do
  len=`cat gene_length.txt | awk '$1~/^'$g'$/' | awk '{print $NF}' | sort -rn | head -n 1`
  awk '$4~/^'$g'$/' gene_exons.bed > $g.bed
  strand=`head -n 1 $g.bed | awk '/-$/' | wc -l`
  if [ $strand -eq 0 ]; then
    # fwd
    ORD="rk4,4"
  else
    # rvs
    ORD="k4,4"
  fi
  bedtools intersect -a $bam -b $g.bed | samtools view -F 256 - | $dir/probe_weight.pl $g 14 $len $strand | awk '{printf $0"	PW:i:"NR"\n"}' | sort -nk1,1 -n$ORD | awk '{printf $0"	TP:i:"NR"\n"}'
  rm $g.bed
done
