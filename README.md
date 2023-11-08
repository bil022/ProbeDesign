## Probe Design
Probe design is the first step of super-resolution FISH experiments to image and quantify DNAs/RNAs. 
Our method is substantially faster and easier to use. We used commonly used 17-nucleotides(nt) counting table to calculate the off-targets. It utilizes a memory-efficient data structure that significantly reduces the memory requirement from 64GB to 2GB. In addition, the method provided multiple parameters and provided new features compared to existing methods.

## Installation

### Dependencies

Linux with minimum 32GB memory is recommended to run the software. Samtools and BWA are required.

* [samtools](https://www.htslib.org/) [1.3.1+]
* [bwa](https://github.com/lh3/bwa) [0.7.12+]

### User installation

Download the source code from Github.

```sh
git clone https://github.com/bil022/ProbeDesign
cd ProbeDesign
```

## Compile source code

To compile from the source code, the library of libboost-regex is required.

```
sudo apt-get update -y
sudo apt-get install -y libboost-regex-dev
```
Compile the source code:

```
make -f bin/Makefile 
```

> g++ -DMAX_REP_SUM=35 -o bin/ProbeDesign -lboost_regex bin/ProbeDesign.cc<br/>
> g++ bin/ProbeDesign.cc -DMAX_REP_SUM=35 -o bin/ProbeDesignStatic -lboost_regex -static -lpthread<br/>
> g++ bin/merfish_pileup.cc -o bin/merfish_pileup

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

## Inputs

* The gene_list.txt include gene of interests

```
cat gene_list.txt 
```

> SLC26A7<br/>
> FMOD<br/>
> SULF1<br/>

* Make sure all gene names are valid, error message will show up for invalid gene names:

```
make -f bin/Makefile id
```
> bin/gene_name.pl gene_list.txt ref/hg38_genes.txt

* Run pipeline:

```
make -f bin/Makefile run
```

The main pipeline for probe design. It will run the following tasks step by step: 

    make -f bin/Makefile gene_list
    make -f bin/Makefile scan_genes
    make -f bin/Makefile gene_qual_design
    make -f bin/Makefile STAR_QUAL
    make -f bin/Makefile filter
In detail, **gene\_list** count the nuclotide frequencies base on the exon information. **scan_genes** count the off-targets for all probe candidates. **gene\_qual\_design** design probes using off-targets information. **STAR\_QUAL** re-map probe sequence to the genome. **filter** identify the top probes.

## Output

The output is **gene_filt48.txt**, the bam and bigbed files are also available for visualization:

```
cat gene_filt48.txt
```
> ...<br/>
> FMOD.14.0.2.104.45      AACTTTTCAGAGAGTGACCACGTCCCTCTG  chr1:203351042:16:30M:AACTTTTCAGAGAGTGACCACGTCCCTCTG<br/>
> SULF1.16.0.2.1.30       TGTCTCTTACGGATTGGAGGCTGAGGGCAG  chr8:69466665:0:30M:CTGCCCTCAGCCTCCAATCCGTAAGAGACA<br/>
> ...<br/>

The **FMOD** and **SULF1** are gene names, **AACTTTTCAGAGAGTGACCACGTCCCTCTG** is the probe sequence, and **chr1:203351042:16** are chromosome and positions, where **0** or **16** present the forward and reverse strand. 

## Esitimated running time
For the installation, most of the tasks require less than 1 min. The building of genome indices (BWS) may take a few hours. For the probe design, one design of 326 genes took about 2 hours. And the job runs in a *batch mode*, it will still need about 1 hour even with one gene. 