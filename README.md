## Probe Design
Probe design is the first step of super-resolution FISH experiments to image and quantify DNAs/RNAs. 
We developed a method with standard 17-nucleotides(nt) off-targets counting table. Our method is substantially faster and easier to use. Notably, it utilizes a memory-efficient data structure that significantly reduces the memory requirement from 64GB to 2GB. In addition, the method provides multiple parameters and new features compared to other methods.

## Installation

### Dependencies

Linux is recommended to run the software. Also, Samtools and BWA is required. 

* [samtools](https://www.htslib.org/)
* [bwa](https://github.com/lh3/bwa)

Also, the index of reference genome (ex. hg38) for bwa should be build if it does not exist.

```
cd ref
bwa index hg38.fa hg38
```

### User installation

Download the source code from Github.

```sh
git clone https://github.com/bil022/ProbeDesigner
cd ProbeDesigner
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

Multiple files (reference genome sequence, repeated sequences, gtf files, indices) are needed for each reference genome. The reference files for human genome are available below, other refrence files are availble upon request:

* [hg38.tgz](http://renlab.sdsc.edu/ProbeDesigner/ref/hg38.tgz/)

## Inputs

Following files are needed as inputs (ex. hg38 for human genome):

* Unzip the required files for reference genomes

```
cd ref
wget http://renlab.sdsc.edu/ProbeDesigner/ref/hg38_ref.tgz
zcat hg38.tgz | tar -xvf -
cd ..
```

* The BWA indices is available in the **ref/hg38**, you can build your own using:

```
cd ref
bwa index hg38.fa hg38
cd ..
```

* Prepare gene_list.txt with gene of interests

```
cat gene_list.txt 
```

> SLC26A7<br/>
> FMOD<br/>
> SULF1<br/>

* Make sure all genes are available, error message will show up for missing gene names:

```
make -f bin/Makefile id
```
> bin/gene_name.pl gene_list.txt ref/hg38_genes.txt

* Run pipeline:

```
make -f bin/Makefile run
```

Main pipeline for probe design, the output is available in **gene_filt48.txt**:

For example:

```
cat gene_filt48.txt
```
> ...<br/>
> FMOD.14.0.2.104.45      AACTTTTCAGAGAGTGACCACGTCCCTCTG  chr1:203351042:16:30M:AACTTTTCAGAGAGTGACCACGTCCCTCTG<br/>
> SULF1.16.0.2.1.30       TGTCTCTTACGGATTGGAGGCTGAGGGCAG  chr8:69466665:0:30M:CTGCCCTCAGCCTCCAATCCGTAAGAGACA<br/>
> ...<br/>

The **FMOD** and **SULF1** are gene names, **AACTTTTCAGAGAGTGACCACGTCCCTCTG** is the probe sequence, and **chr1:203351042:16** are chromosome and positions, where **0** or **16** present the forward and reverse strand. 