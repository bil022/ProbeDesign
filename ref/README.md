## Download files for reference genomes

Multiple files of reference genome sequence, gene annotation files (.gtf) and counting table files (.idx) are needed for the reference genome. The reference files for human genome are available below, other refrence files are availble upon request:

* [hg38_ref.tgz](http://renlab.sdsc.edu/ProbeDesigner/ref/hg38_ref.tgz/)

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

