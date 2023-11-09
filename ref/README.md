## Download files for reference genomes

Multiple files of reference genome sequence, gene annotation files (.gtf) and counting table files (.idx) are needed for the reference genome. The reference files for human genome are available below, other refrence files are availble upon request:

* [hg38_ref.tgz](http://renlab.sdsc.edu/ProbeDesigner/ref/hg38_ref.tgz)

Download and unzip hg38 reference genome:

```
cd ref
wget http://renlab.sdsc.edu/ProbeDesigner/ref/hg38_ref.tgz
zcat hg38.tgz | tar -xvf -
cd ..
```

* The BWA indices should be available in the **ref/hg38**, to build the indices:

```
cd ref
bwa index hg38.fa hg38
cd ..
```

* To build your own reference genome files:

Following is the script to build hg38 reference genome files, you can use it to build other reference genome:

        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
        zcat hg38.ncbiRefSeq.gtf.gz | grep -v ^chrM | awk '$1!~/_/' | gffread -W -w hg38_transcriptW.fa -g hg38.fa -
        ../bin/ProbeDesign build_rna hg38_transcriptW.fa hg38_transcriptW.idx 2>&1 | gzip > hg38_build_rna.log.gz
        zcat hg38.ncbiRefSeq.gtf.gz | grep -v ^chrM | grep -Pe "\t(exon|3UTR)\t" | awk '$1!~/_/' | sort -k10,10 > hg38_exons.gtf
        ../bin/gtf2bed.pl hg38_exons.gtf | ../bin/merfish_pileup pileup | ../bin/merfish_pileup pileup_fa hg38.fa > hg38_exons.fa 2> hg38.chrom.sizes 
        ../bin/gtf2bed.pl hg38_exons.gtf | ../bin/merfish_pileup pileup | sed 's/:/ /g' | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$3,$6,$2}' | sort -k1,1 -k2,2n > hg38_exons.bed
        awk '/gene=/{print $2}' hg38_transcriptW.fa | sort -u | sed 's/=/ /' | awk '{print $NF}' > hg38_genes.txt
