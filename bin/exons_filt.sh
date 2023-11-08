#!/bin/bash

dir=`dirname $0`
bam=gene_design_uniq.bam
samtools view -H $bam
$dir/exons.pl $dir/../ref/hg38_exons.bed > gene_exons.bed
for g in `awk '{print $4}' gene_exons.bed | uniq`; do
  len=`cat gene_length.txt | awk '$1~/^'$g'$/' | awk '{print $NF}' | sort -rn | head -n 1`
  awk '$4~/^'$g'$/' gene_exons.bed > $g.bed
  bedtools intersect -a $bam -b $g.bed | samtools view -F 256 - | $dir/probe_weight.pl $g 14 $len | head -n 48
  rm $g.bed
done
