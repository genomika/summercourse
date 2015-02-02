[WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Alinhamento__
% _(updated 29-01-2014)_

<!-- COMMON LINKS HERE -->

[BWA]: http://bio-bwa.sourceforge.net/ "BWA"
[HPG Aligner]: https://github.com/opencb/hpg-aligner/wiki/ "HPG Aligner"
[Bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml "Bowtie2"
[TopHat]: http://ccb.jhu.edu/software/tophat/index.shtml "TopHat2"
[STAR]: https://code.google.com/p/rna-star/ "STAR"
[MapSplice2]: http://www.netlab.uky.edu/p/bioinfo/MapSplice2 "MapSplice2"
[SAMTools]: http://samtools.sourceforge.net/ "SAMtools"
[dwgsim]: http://sourceforge.net/apps/mediawiki/dnaa/index.php?title=Whole_Genome_Simulation "dwgsim"
[BEERS]: http://www.cbil.upenn.edu/BEERS/ "BEERS"
[Ensembl]: http://www.ensembl.org/index.html "Ensembl"
[Picard]: http://broadinstitute.github.io/picard/ "Picard"

# Preliminaries

In this hands-on will learn how to align DNAseq data with most widely used software today. Building a whole genome index requires a lot of RAM memory and almost one hour in a typical workstation, for this reason **in this tutorial we will work with targeted panel of genes data** to speed up the exercises. The same steps would be done for a whole genome alignment.

Once raw sequence files are generated (in FASTQ format) and quality-checked, the next step in most NGS pipelines is mapping to a reference genome. For individual sequences of interest, it is common to use a tool like BLAST to identify genes or species of origin. However, a typical example NGS dataset may have tens to hundreds of millions of reads, which BLAST and similar tools are not designed to handle.
Thus, a large set of computational tools have been developed to quickly, and with sufficient (but not absolute) accuracy align each read to its best location, if any, in a reference. Even though many mapping tools exist, a few individual programs have a dominant "market share" of the NGS world. These programs vary widely in their design, inputs, outputs, and applications. In this section, we will primarily focus on two of the most versatile mappers: BWA and Bowtie2 (only comments).

### NGS aligners discussed:

- [BWA] : BWA is a software package for mapping **DNA** low-divergent sequences against a large reference genome, such as the human genome.
- [Bowtie2] : *Bowtie 2* is an ultrafast and memory-efficient tool for aligning **DNA** sequencing reads to long reference sequences (only comments).

### Other software used in this hands-on:
- [SAMTools] : SAM Tools **provide various utilities** for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [Picard] : Piccard Tools **provide various utilities** for handling alignments in the SAM and BAM format, including sorting, merging, indexing and generating alignments in a per-position format.

### File formats explored:

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf): Sequence alignment format, plain text.
- [BAM](http://www.broadinstitute.org/igv/bam): Binary and compressed version of SAM


### Data used in this practical

You have already worked with two human data samples, which we will continue to use here.  The paired end data should already be located at:

    $WORK/archive/original/2014_05.core_ngs

- Sample_Yeast_L005_R1.cat.fastq.gz	Paired-end Illumina, First of pair, FASTQ	Yeast ChIP-seq
- Sample_Yeast_L005_R2.cat.fastq.gz	Paired-end Illumina, Second of pair, FASTQ	Yeast ChIP-seq

First copy the two human datasets to your $SCRATCH/core_ngs/fastq_prep directory.

    cd $SCRATCH/core_ngs/fastq_prep
    cp $CLASSDIR/human_stuff/*rnaseq.fastq.gz .

Create a $SCRATCH/core_ngs/align directory and make a link to the fastq_prep directory.

    mkdir -p $SCRATCH/core_ngs/align
    cd $SCRATCH/core_ngs/align
    ln -s -f ../fastq_prep fq
    ls -l
    ls fq

#### Reference Genomes

Before we get to alignment, we need a genome to align to.  We will use the human genome (hg19) here.

    hg19	Human	3.1 Gbp	25 (really 93)	UCSC	UCSC GoldenPath

Searching genomes is hard work and takes a long time if done on an un-indexed, linear genomic sequence.  So aligners require that references first be indexed for quick access  The aligners we are using each require a different index, but use the same method (the Burrows-Wheeler Transform) to get the job done. This requires taking a FASTA file as input, with each chromosome (or contig) as a separate entry, and producing some aligner-specific set of files as output. Those index files are then used by the aligner when performing the sequence alignment. 

**hg19** is way too big for us to index here, so we're not going to do it. Instead, all hg19 index files are located at:

    /scratch/01063/abattenh/ref_genome/bwa/bwtsw/hg19


First stage the yeast and mirbase reference FASTA files in your work archive area in a directory called references.

    mkdir -p $WORK/archive/references/fasta
    cp $CLASSDIR/references/*.fa $WORK/archive/references/fasta/

 With that, we're ready to get started on the first exercise.




Create a ```data``` folder in your **working directory** and download the **reference genome sequence** to be used (human chromosome 21) and *simulated datasets* from **Dropbox** [data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a).
For the rest of this tutorial the **working directory** will be **cambridge_mda14** and all the **paths** will be relative to that working directory:
    
    cd cambridge_mda14
    mkdir data

##### Download reference genome from [Ensembl]

Working with NGS data requires a high-end workstations and time for building the reference genome indexes and alignment. During this tutorial we will work only with chromosome 21 to speed up the runtimes. You can download it from **Dropbox** [data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a) or from the *Download* link at the top of [Ensembl] website and then to *Download data via FTP*, you get it in only one step by going to:

[http://www.ensembl.org/info/data/ftp/index.html](http://www.ensembl.org/info/data/ftp/index.html)

You should see a species table with a Human (*Homo sapiens*) row and a *DNA (FASTA)* column or click at [ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/), download the chromosome 21 file (*Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz*) and move it from your browser download folder to your ```data``` folder:

    mv Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz path_to_local_data

**NOTE:** For working with the whole reference genome the file to be downloaded is *Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz*


##### Copy simulated datasets

For this hands-on we are going to use small DNA and RNA-seq datasets simulated from chromosome 21. Data has been already simulated using _dwgsim_ software from SAMtools for DNA and _BEERS_ for RNA-seq. You can copy from the shared resources from **Dropbox** [data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a) into your ``data`` directory for this practical session. Preparing the data directory:

    cp path_to_course_materials/alignment/* your_local_data/

Notice that the name of the folders and files describe the dataset, ie. ```dna_chr21_100_hq``` stands for: _DNA_ type of data from _chromosome 21_ with _100_nt read lengths of _high_ quality. Where _hq_ quality means 0.1% mutations and _lq_ quality 1% mutations. Take a few minutes to understand the different files.

**NOTE:** If you want to learn how to simulate DNA and RNA-seq for other conditions go down to the end of this tutorial.

##### Real datasets

For those with access to high-end nodes clusters you can index and simulated whole genome datasets or download real datasets from this sources:
- [1000genomes project](http://www.1000genomes.org/)
- [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/)
- [Sequence Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra)


### Installing SAMtools (Optional, already installed)

Check that it is not installed by executing

    samtools

A list of commands should be printed. If not then proceed with the installation.

Download [SAMtools] from *SF Download Page* link and move to the working directory, then uncompress it.

    mv samtools-0.1.19.tar.bz2 working_directory
    cd working_directory
    tar -jxvf samtools-0.1.19.tar.bz2 
    cd samtools-0.1.19
    make

Check that is correct by executing it with no arguments, the different commands available should be printed. You can also copy it to your ```bin``` folder in your home directory, if bin folder exist, to make it available to the PATH:

    samtools
    cp samtools ~/bin


# Exercise 1: NGS Genomic DNA aligment

In this exercise we'll learn how to download, install, build the reference genome index and align in single-end and paired-end mode with the two most widely DNA aligners: *BWA* and *Bowtie2*. But first, create an ```aligners``` folder to store the software, and an ```alignments``` folder to store the results, create those folders in your *working directory* next to ```data```, you can create both folders by executing:

    mkdir aligners
    mkdir alignments

Now go to ```aligners``` and  ```alignments``` folders and create subfolders for *bwa* and *bowtie* to store the indexes and alignments results:

    cd aligners
    mkdir bwa hpg-aligner bowtie

and

    cd alignments
    mkdir bwa hpg-aligner bowtie
    
**NOTE:** Now your working directory must contain 3 folders: data (with the reference genome of chrom. 21 and simulated datasets), aligners and alignments. Your working directory should be similar to this (notice that aligners have not been downloaded):

```
.
├── aligners
│   ├── bowtie
│   ├── bwa
│   ├── hpg-aligner
├── alignments
│   ├── bowtie
│   ├── bwa
│   ├── hpg-aligner
├── data
│   ├── dna_chr21_100_hq_read1.fastq
│   ├── dna_chr21_100_hq_read2.fastq
│   ├── dna_chr21_100_lq_read1.fastq
│   ├── dna_chr21_100_lq_read2.fastq
│   ├── Homo_sapiens_cDNAs_chr21.fa
│   ├── Homo_sapiens.GRCh37.75.dna.chromosome.21.fa
│   ├── rna_chr21_100_hq_read1.fastq
│   └── rna_chr21_100_hq_read2.fastq
│   ├── rna_chr21_150_lq_read1.fastq
│   └── rna_chr21_150_lq_read2.fastq


```


### BWA
[BWA] is probably the most used aligner for DNA. AS the documentation states it consists of three different algorithms: *BWA*, *BWA-SW* and *BWA-MEM*. The first algorithm, which is the oldest, is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA for 70-100bp Illumina reads.

All these three algorithms come in the same binary so only one download and installation is needed.

##### Download and install (Optional, already installed)

First check that bwa is not currently installed by executing:

    bwa

A list of commands will be printed if already installed. If not you can continue with the installation.

You can click on ```SF download page``` link in the [BWA] page or click directly to:

[http://sourceforge.net/projects/bio-bwa/files](http://sourceforge.net/projects/bio-bwa/files)

Click in the last version of BWA and wait for a few seconds, as the time of this tutorial last version is **bwa-0.7.10.tar.bz2**, the download will start. When downloaded go to your browser download folder and move it to aligners folder, uncompress it and compile it:

    mv bwa-0.7.10.tar.bz2 working_directory/aligners/bwa
    tar -jxvf bwa-0.7.10.tar.bz2
    cd bwa-0.7.10
    make
    cp bwa ~/bin

You can check that everything is allright by executing:

    bwa

Some information about the software and commands should be listed.


##### Build the index

Create a folder inside ```aligners/bwa``` folder called ```index``` to store the BWA index and copy the reference genome into it:
    
    mkdir index
    cp  data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa  aligners/bwa/index/   (this path can be different!)
    
Now you can create the index by executing:

    bwa index aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa

Some files will be created in the ```index``` folder, those files constitute the index that BWA uses.

**NOTE:** The index must created only once, it will be used for all the different alignments with BWA.


##### Aligning with new BWA-MEM in both single-end (SE) and paired-end (PE) modes

BWA-MEM is the recommended algorithm to use now. You can check the options by executing:

    bwa mem

To align **SE** with BWA-MEM execute:

    bwa mem -t 4 -R "@RG\tID:foo\tSM:bar\tPL:Illumina\tPU:unit1\tLB:lib1" aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa data/dna_chr21_100_hq_read1.fastq > alignments/bwa/dna_chr21_100_hq_se.sam

Now you can use SAMtools to create the BAM file from the *alignment/bwa* folder:
    cd alignments/bwa
    samtools view -S -b dna_chr21_100_hq_se.sam -o dna_chr21_100_hq_se.bam


To align **PE** with BWA-MEM just execute the same command line with the two FASTQ files:

    bwa mem -t 4 -R "@RG\tID:foo\tSM:bar\tPL:Illumina\tPU:unit1\tLB:lib1" aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa data/dna_chr21_100_hq_read1.fastq data/dna_chr21_100_hq_read2.fastq > alignments/bwa/dna_chr21_100_hq_pe.sam
    
Now you can use SAMtools to create the BAM file from the *alignment/bwa* folder:

    cd alignments/bwa
    samtools view -S -b dna_chr21_100_hq_pe.sam -o dna_chr21_100_hq_pe.bam


Now you can do the same for the **low** quality datasets.


##### Aligning with old BWA algorithm, two command lines: ALN and SAMSE/SAMPE in SE and PE modes (Optional exercise)

Now we are going to align SE and PE the **high** quality dataset. Single-end alignment with BWA requires 2 executions. The first uses ```aln``` command and takes the ```fastq``` file and creates a ```sai``` file; the second execution uses ```samse``` and the ```sai``` file and create the ```sam``` file. Results are stored in ```alignments``` folder:

    bwa aln aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 data/dna_chr21_100_hq_read1.fastq -f alignments/bwa/dna_chr21_100_hq_se.sai
    bwa samse aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa alignments/bwa/dna_chr21_100_hq_se.sai data/dna_chr21_100_hq_read1.fastq -f alignments/bwa/dna_chr21_100_hq_se.sam


For paired-end alignments with BWA 3 executions are needed: 2 for ```aln``` command and 1 for ```sampe``` command:

    bwa aln aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 data/dna_chr21_100_hq_read1.fastq -f alignments/bwa/dna_chr21_100_hq_pe1.sai
    bwa aln aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 data/dna_chr21_100_hq_read2.fastq -f alignments/bwa/dna_chr21_100_hq_pe2.sai
    bwa sampe aligners/bwa/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa alignments/bwa/dna_chr21_100_hq_pe1.sai alignments/bwa/dna_chr21_100_hq_pe2.sai data/dna_chr21_100_hq_read1.fastq data/dna_chr21_100_hq_read2.fastq -f alignments/bwa/dna_chr21_100_hq_pe.sam

Now you can use SAMtools to create the BAM file from the *alignment/bwa* folder:

    cd alignments/bwa
    samtools view -S -b dna_chr21_100_hq_se.sam -o dna_chr21_100_hq_se.bam
    samtools view -S -b dna_chr21_100_hq_pe.sam -o dna_chr21_100_hq_pe.bam

Now you can do the same for the **low** quality datasets.


### HPG Aligner
[HPG Aligner] is an ultrafast and high sensitivity tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of any size. HPG Aligner 2.0 indexes the genome using a Suffix Arrays. HPG Aligner 2 supports gapped, local, and paired-end alignment modes, together with INDEL realignment and recalibration.

##### Build the index

Create a folder inside HPG Aligner program called ```index``` to store the HPG Aligner index and copy the reference genome into it:
    
    cd hpg-aligner   (_inside aligners folder_)
    mkdir index
    cp data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa aligners/bwa/index/
    
Now you can create the index by executing:

    hpg-aligner build-sa-index -g aligners/hpg-aligner/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -i aligners/hpg-aligner/index/

Some files will be created in the ```index``` folder, those files constitute the index that HPG Aligner uses.

**NOTE:** The index must created only once, it will be used for all the different alignments with HPG Aligner.

### Bowtie2

[Bowtie2] as documentation states is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to few 100s. Bowtie 2 indexes the genome with an FM Index to keep its memory footprint small: for the human genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.

##### Download and install (Optional, already installed)

First check that bwa is not currently installed by executing:

    bowtie2

A list of commands will be printed if already installed. If not you can continue with the installation.

From [Bowtie2] go to ```Latest Release``` and download the program or go directly to:

[http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/)

Click in the Linux version of Bowtie2 and wait for a few seconds, as the time of this tutorial last version is **bowtie2-2.2.3-linux-x86_64.zip**, the download will start. When downloaded go to your browser download folder and move it to aligners folder and uncompress it. No need to compile if you downloaded the Linux version:

    mv bowtie2-2.2.3-linux-x86_64.zip working_directory/aligners/bowtie
    unzip bowtie2-2.2.3-linux-x86_64.zip
    cd bowtie2-2.2.3

You can check that everything is allright by executing:

    bowtie2

Big information about the software and commands should be listed.


##### Build the index

Create a folder inside Bowtie2 program called ```index``` to store the Bowtie2 index and copy the reference genome into it:
    
    cd bowtie   (_inside aligners folder_)
    mkdir index
    cp data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa aligners/bowtie/index/
    
Now you can create the index by executing:

    bowtie2-build aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa

Some files will be created in the ```index``` folder, those files constitute the index that Bowtie2 uses.

**NOTE:** The index must created only once, it will be used for all the different alignments with Bowtie2.


##### Aligning in SE and PE modes

Mapping **SE** with Bowtie2 requires only 1 execution, for aligning the **high** in SE mode execute:

    bowtie2 -q -p 4 -x aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -U data/dna_chr21_100_hq_read1.fastq -S alignments/bowtie/dna_chr21_100_hq_se.sam

And create the BAM file using SAMtools;

    cd alignments/bowtie
    samtools view -S -b dna_chr21_100_hq_se.sam -o dna_chr21_100_hq_se.bam


Mapping in **PE** also requires only one execution:

    bowtie2 -q -p 4 -x aligners/bowtie/index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -1 data/dna_chr21_100_hq_read1.fastq -2 data/dna_chr21_100_hq_read2.fastq -S alignments/bowtie/dna_chr21_100_hq_pe.sam

And create the BAM file using SAMtools;

    cd alignments/bowtie
    samtools view -S -b dna_chr21_100_hq_pe.sam -o dna_chr21_100_hq_pe.bam
    
Repeat the same steps for the **low** quality dataset.


### More exercises

- Try to simulate datasets with longer reads and more mutations to study which aligner behaves better
- Test the aligner sensitivity to INDELS
- Try BWA-MEM algorithm and compare sensitivity. The same index is valid, only one execution for the SAM file ```./bwa mem index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read1.fastq```
