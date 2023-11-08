ref = hg38

all:	
	g++ -DMAX_REP_SUM=35 -o bin/ProbeDesign -lboost_regex bin/ProbeDesign.cc
	g++ bin/ProbeDesign.cc -DMAX_REP_SUM=35 -o bin/ProbeDesignStatic -lboost_regex -static -lpthread
	g++ bin/merfish_pileup.cc -o bin/merfish_pileup
id:
	bin/gene_name.pl gene_list.txt ref/hg38_genes.txt
run:
	$(MAKE) -f bin/Makefile gene_list
	$(MAKE) -f bin/Makefile scan_genes
	$(MAKE) -f bin/Makefile gene_qual_design
	$(MAKE) -f bin/Makefile STAR_QUAL
	$(MAKE) -f bin/Makefile filter
gene_list:
	cat ref/$(ref)_transcriptW.fa | bin/transcriptW.pl `cat gene_list.txt` 2> gene_length.txt > gene_list.fastq
scan_genes:
	bin/ProbeDesignStatic scan ref/$(ref)_transcriptW.idx ref/$(ref)_transcriptW_rep.idx gene_list.fastq gene_list.cnts
gene_qual_design:
	bin/ProbeDesignStatic design 14 gene_list.cnts 2> gene_qual_pool.txt > gene_qual_design.txt
	bin/gene_id.pl gene_qual_design.txt | sort -u > gene_qual_design_uniq.txt
	awk '{print $$1}' gene_qual_design_uniq.txt | uniq -c | awk '{print $$2"\t"$$1}' > gene_qual_design_cnts.txt
STAR_QUAL:
	awk '{print ">"$$1"."$$3"."$$4"."$$5"\n"$$2;}' gene_qual_design_uniq.txt > gene_qual_design_uniq.fa
	STAR --genomeDir ref/$(ref)/ --readFilesIn gene_qual_design_uniq.fa
	samtools view -Sb Aligned.out.sam | samtools sort - > gene_qual_design_uniq.bam
	samtools index gene_qual_design_uniq.bam
filter:
	bin/exons.pl ref/$(ref)_exons.bed > gene_exons.bed
	bin/exons_all.sh | samtools view -Sb - | samtools sort - > gene_all.bam
	samtools index gene_all.bam 
	samtools view -h gene_all.bam | bin/SamTagFilt.pl TP | samtools view -Sb - > gene_3prime.bam
	samtools index gene_3prime.bam
	samtools view gene_all.bam|sed 's/\./ /' | awk '{print $$1}' | sort | uniq -c > gene.all_cnt.txt
	#bin/exons_filt.sh | samtools view -Sb - | samtools sort - > gene_filt48.bam
	samtools view -h gene_all.bam | bin/limits.pl | samtools view -Sb - > gene_filt48.bam
	samtools index gene_filt48.bam 
	samtools view gene_filt48.bam|sed 's/\./ /' |awk '{print $$1}' | uniq -c > gene_filt48.cnt.txt
	bamToBed -split -i gene_filt48.bam | sort -k1,1 -k2,2n > gene_filt48.bed
	bedToBigBed gene_filt48.bed ref/$(ref).chrom.sizes gene_filt48.bb
	samtools view gene_filt48.bam|awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$3,$$4,$$6,$$10}' | bin/bam2txt.pl > gene_filt48.txt
mail:
	tar -cvzf MERFISH_src.tgz Makefile *.sh *.pl Readme.txt 
	pwd | mail -s MERFISH_SRC -a MERFISH_src.tgz bil022@ucsd.edu
tomato:
	bin/ProbeDesignStatic scan ref/${ref}_transcriptW.idx ref/${ref}_transcriptW_rep.idx tdTomato.fa tdTomato.cnts
	bin/ProbeDesignStatic design 14 tdTomato.cnts 2> tdTomato_pool.txt > tdTomato.txt

