[WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Alinhamento__
% _(updated 29-01-2014)_

<!-- COMMON LINKS HERE -->

[SAMTools]: http://samtools.sourceforge.net/ "SAMtools"
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

    workshop@172.16.225.11:/home/workshop/ngscourse/data/sample2/

The files you will be working with:

- 510-7-BRCA_S8_L001_R1.fastq.gz Paired-end Illumina, First of pair, FASTQ Human Sample-seq
- 510-7-BRCA_S8_L001_R2.fastq.gz Paired-end Illumina, Second of pair, FASTQ Human Sample-seq

First copy the two human dataset to your $HOME/core_ngs/fastq_prep directory.

    cd $HOME/core_ngs/fastq_prep
    cp $HOME/work/data/sample1/sample1/*.fastq.gz .
    cp $HOME/work/data/sample2/*.fastq.gz .

Create a $HOME/core_ngs/align directory and make a link to the fastq_prep directory.

    mkdir -p $HOME/core_ngs/align
    cd $HOME/core_ngs/align
    ln -s -f ../fastq_prep fq
    ls -l
    ls fq

#### Reference Genomes

Before we get to alignment, we need a genome to align to.  We will use the human genome (hg19) here.

    hg19	Human	3.1 Gbp	25 (really 93)	UCSC	UCSC GoldenPath

Searching genomes is hard work and takes a long time if done on an un-indexed, linear genomic sequence.  So aligners require that references first be indexed for quick access  The aligners we are using each require a different index, but use the same method (the Burrows-Wheeler Transform) to get the job done. This requires taking a FASTA file as input, with each chromosome (or contig) as a separate entry, and producing some aligner-specific set of files as output. Those index files are then used by the aligner when performing the sequence alignment. 

**hg19** is way too big for us to index here, so we're not going to do it. Instead, all hg19 index files are located at:

    $HOME/data


 With that, we're ready to get started on the first exercise.

Exercise #1: BWA aln 

We will perform a global alignment of the paired-end Human sequences using bwa. This workflow generally has the following steps:


- Trim the FASTQ sequences down to 50 with fastx_clipper
   - this removes most of any 5' adapter contamination without the fuss of specific adapter trimming w/cutadapt
- Prepare the reference index for bwa (one time) using bwa index
- Perform a global bwa alignment on the R1 reads (bwa aln) producing a BWA-specific binary .sai intermediate file
- Perform a global bwa alignment on the R2 reads (bwa aln) producing a BWA-specific binary .sai intermediate file
- Perform pairing of the separately aligned reads and report the alignments in SAM format using bwa sampe
- Convert the SAM file to a BAM file (samtools view)
- Sort the BAM file by genomic location (samtools sort or use Picard tool to generate the bam)
- Index the BAM file (samtools index)
- Gather simple alignment statistics (samtools flagstat and samtools idxstat)

We're going to skip the trimming step for now and see how it goes. We'll perform steps 2 - 5 now and leave samtools for the next course section, since those steps (6 - 10) are common to nearly all post-alignment workflows.

### Introducing BWA

Like other tools you've worked with so far, you first need to load bwa using the module system.  Go ahead and do that now, and then enter bwa with no arguments to view the top-level help page (many NGS tools will provide some help when called with no arguments).

    bwa
    

    Program: bwa (alignment via Burrows-Wheeler transformation)
    Version: 0.7.7-r441

    Contact: Heng Li <lh3@sanger.ac.uk>
    Usage:   bwa <command> [options]
    
    Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
         pemerge       merge overlapping paired ends (EXPERIMENTAL)
         aln           gapped/ungapped alignment
         samse         generate alignment (single ended)
         sampe         generate alignment (paired ended)
         bwasw         BWA-SW for long queries
         fa2pac        convert FASTA to PAC format
         pac2bwt       generate BWT from PAC
         pac2bwtgen    alternative algorithm for generating BWT
         bwtupdate     update .bwt to the new format
         bwt2sa        generate SA from BWT and Occ
    Note: To use BWA, you need to first index the genome with `bwa index'.
      There are three alignment algorithms in BWA: `mem', `bwasw', and
      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
      first. Please `man ./bwa.1' for the manual. 


As you can see, bwa offers a number of sub-commands one can use with to do different things.

#### Building the BWA  hg19 index

We're going to index the genome with the index command. To learn what this sub-command needs in the way of options and arguments, enter bwa index with no arguments.

    Usage:   bwa index [-a bwtsw|is] [-c] <in.fasta>
    Options: -a STR    BWT construction algorithm: bwtsw or is [auto]
          -p STR    prefix of the index [same as fasta name]
          -6        index files named as <in.fasta>.64.* instead of <in.fasta>.*
    Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
         `-a div' do not work not for long genomes. Please choose `-a'
         according to the length of the genome.


Here, we only need to specify two things:
-   the name of the FASTA file  
-   whether to use the  bwtsw or is algorithm for indexing

Since the genome is relative large (~ 3.1 Gb) we will specify bwtsw as the indexing option, and the name of the FASTA file is ucsc.hg19.fasta.

Importantly, the output of this command is a group of files that are all required together as the index.

Now execute the bwa index command.

    bwa index -a bwtsw ucsc.hg19.fasta

Since the human genome is quite large, this should  take long to execute.  We already provide for you the genome indexed. When it is complete you should see a set of index files like this:

    ucsc.hg19.fasta
    ucsc.hg19.amb
    ucsc.hg19.ann
    ucsc.hg19.bwt
    ucsc.hg19.pac
    ucsc.hg19.sa


A common question is what contigs are in a given FASTA file. You'll usually want to know this before you start the alignment so that you're familiar with the contig naming convention – and to verify that it's the one you expect.

We saw that a FASTA consists of a number of contig entries, each one starting with a name line of the form below, followed by many lines of bases.

    >contigName

How do we dig out just the lines that have the contig names and ignore all the sequences? Well, the contig name lines all follow the pattern above, and since the > character is not a valid base, it will never appear on a sequence line.
We've discovered a pattern (also known as a regular expression) to use in searching, and the command line tool that does regular expression matching is **grep**.

Regular expressions are so powerful that nearly every modern computer language includes a "regex" module of some sort. There are many online tutorials for regular expressions (and a few different flavors of them). But the most common is the [Perl style](http://perldoc.perl.org/perlretut.html). We're only going to use the most simple of regular expressions here, but learning more about them will pay handsome dividends for you in the future.

Here's how to execute grep to list contig names in a FASTA file.

    grep -P '^>' ucsc.hg19.fasta | more

**Notes:**
- The -P option tells grep to use Perl-style regular expression patterns. 
- This makes including special characters like Tab ( \t ), Carriage Return ( \r ) or Linefeed ( \n ) much easier that the default Posix paterns.
- While it is not really required here, it generally doesn't hurt to include this option.
- '^>' is the regular expression describing the pattern we're looking for (described below)
- ucsc.hg19.fasta is the file to search. Lines with text that match our pattern will be written to standard output; non matching lines will be omitted.
- We pipe to more just in case there are a lot of contig names.

Now down to the nuts and bolts of our pattern, '^>'

First, the single quotes around the pattern – they are only a signal for the bash shell. As part of its friendly command line parsing and evaluation, the shell will often look for special characters on the command line that mean something to it (for example, the $ in front of an environment variable name, like in $HOME). Well, regular expressions treat the $ specially too – but in a completely different way! Those single quotes tell the shell "don't look inside here for special characters – treat this as a literal string and pass it to the program". The shell will obey, will strip the single quotes off the string, and will pass the actual pattern, ^>, to the grep program. (Aside: We've see that the shell does look inside double quotes ( " ) for certain special signals, such as looking for environment variable names to evaluate.)

So what does ^> mean to grep? Well, from our contig name format description above we see that contig name lines always start with a > character, so > is a literal for grep to use in its pattern match.

We might be able to get away with just using this literal alone as our regex, specifying '>' as the command line argument. But for grep, the more specific the pattern, the better. So we constrain where the > can appear on the line. The special carat ( ^ ) character represents "beginning of line". So ^> means "beginning of a line followed by a > character, followed by anything. (Aside: the dollar sign ( $ ) character represents "end of line" in a regex. There are many other special characters, including period ( . ), question mark ( ? ), pipe ( | ), parentheses ( ( ) ), and brackets ( [ ] ), to name the most common.)


#### Exercise: How many contigs are there in the human genome reference?

    grep -P -c '^>' ucsc.hg19.fasta


### Performing the bwa alignment

Now, we're ready to execute the actual alignment, with the goal of initially producing a SAM file from the input FASTQ files and reference. First go to the align directory, and link to the human reference directory (this will make our commands more readable).

    cd $HOME/core_ngs/align
    ln -s /mnt/data/
    ls data/

As our workflow indicated, we first use bwa aln on the R1 and R2 FASTQs, producing a BWA-specific .sai intermediate binary files. Since these alignments are completely independent, we can execute them in parallel in a batch job.

There are lots of options, but here is a summary of the most important ones. BWA, is a lot more complex than the options let on. If you look at the BWA manual on the web for the aln sub-command, you'll see numerous options that can increase the alignment rate (as well as decrease it), and all sorts of other things. 

-l	Controls the length of the seed (default = 32)
-k	Controls the number of mismatches allowable in the seed of each alignment (default = 2)
-n	Controls the number of mismatches (or fraction of bases in a given alignment that can be mismatches) in the entire alignment (including the seed) (default = 0.04)
-t	Controls the number of threads

The rest of the options control the details of how much a mismatch or gap is penalized, limits on the number of acceptable hits per read, and so on.  Much more information can be accessed at the [BWA manual page](http://bio-bwa.sourceforge.net/bwa.shtml).

For a simple alignment like this, we can just go with the default alignment parameters, with one exception. At the server, we want to optimize our alignment speed by allocating more than one thread (-t) to the alignment. We want to run 2 tasks, and will use a minimum of one 2-core node. So we can assign 2 cores to each alignment by specifying -t 2.

Also note that bwa writes its (binary) output to standard output by default, so we need to redirect that to a .sai file. And of course we need to redirect standard error to a log file, one per file. The '&' means run in background. It is useful for long jobs to run in background.

Run the following lines:

    bwa aln -t 2 data/ucsc.hg19 fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R1_001.fastq.gz > BRCA_R1.sai 2> aln.brca_R1.log &
    bwa aln -t 2 data/ucsc.hg19 fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R2_001.fastq.gz > BRCA_R2.sai 2> aln.brca_R2.log &

Since you have directed standard error to log files, you can use a neat trick to monitor the progress of the alignment: tail -f. The -f means "follow" the tail, so new lines at the end of the file are displayed as they are added to the file.

    # Use Ctrl-c to stop the output any time
    tail -f aln.brca_R1.log

When it's done you should see two .sai files. Next we use the bwa sampe command to pair the reads and output SAM format data. For this command you provide the same reference prefix as for bwa aln, along with the two .sai files and the two original FASTQ files.

Again bwa writes its output to standard output, so redirect that to a .sam file. (Note that bwa sampe is "single threaded" – it does not have an option to use more than one processor for its work.) We'll just execute this at the command line.

    bwa sampe data/ucsc.hg19 BRCA_R1.sai BRCA_R2.sai fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R1_001.fastq.gz  fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R2_001.fastq.gz > brca_pairedend.sam

You did it!  You should now have a SAM file that contains the alignments. It's just a text file, so take a look with head, more, less, tail, or whatever you feel like. In the next section, with samtools, you'll learn some additional ways to analyze the data once you create a BAM file.

#### Exercise: What kind of information is in the first lines of the SAM file?

The SAM or BAM has a number of header lines, which all start with an at sign ( @ ). The @SQ lines describe each contig and its length. There is also a @PG  line that describes the way the bwa sampe was performed.

#### Exercise: How many alignment records (not header records) are in the SAM file?

This looks for the pattern  '^M00538' which is the start of every read name (which starts every alignment record).
Remember -c says just count the records, don't display them.

    grep -P -c '^M00538' brca_pairedend.sam

Or use the -v (invert) option to tell grep to print all lines that don't match a particular pattern, here the header lines starting with @.

     grep -P -v -c '^@' brca_pairedend.sam


#### Exercise: How many sequences were in the R1 and R2 FASTQ files combined?
#### Do both R1 and R2 reads have separate alignment records?
#### Does the SAM file contain both aligned and un-aligned reads?
#### What is the order of the alignment records in this SAM file?

**Answers**

- 128552 + 128552 sequences
- yes, it must, because there were 128552 + 128552 R1+R2 reads and an equal number of alignment records
- yes, it must, because there were 128552 + 128552 R1+R2 reads and an equal number of alignment records
- the names occur in the exact same order as they did in the FASTQ, except that they come in pairs
   -   the R1 read comes first, then its corresponding R2
   -   this ordering is called read name ordering


### Using cut to isolate fields

Suppose you wanted to look only at field 3 (contig name) values in the SAM file. You can do this with the handy cut command. Below is a simple example where you're asking cut to display the 3rd of the last 10 alignments.

    tail brca_pairedend.sam | cut -f 3

By default cut assumes the field delimiter is Tab, which is the delimiter used in the majority of NGS file formats. You can, of course, specify a different delimiter with the -d option.

You can also specify a range of fields, and mix adjacent and non-adjacent fields. This displays fields 2 through 6, field 9, and all fields starting with 12 (SAM tag fields).

    tail brca_pairedend.sam | cut -f 2-6,9,12-

You may have noticed that some alignment records contain contig names (e.g. chrV) in field 3 while others contain an asterisk ( * ). Usually the * means the record didn't align. (This isn't always true – later you'll see how to properly distinguish between mapped and unmapped reads using samtools.) We're going to use this heuristic along with cut to see about how many records represent aligned sequences.

First we need to make sure that we don't look at fields in the SAM header lines. We're going to end up with a series of pipe operations, and the best way to make sure you're on track is to enter them one at a time piping to head:

     # the ^M00538 pattern matches lines starting with HWI (the start of all read names in column 1)
    grep -P '^M00538' brca_pairedend.sam | head

Ok, it looks like we're seeing only alignment records. Now let's pull out only field 3 using cut:

    grep -P '^M00538' brca_pairedend.sam | cut -f 3 | head

Cool, we're only seeing the contig name info now. Next we use grep again, piping it our contig info and using the -v (invert) switch to say print lines that don't match the pattern:

    grep -P '^M00538' brca_pairedend.sam | cut -f 3 | grep -v '*' | head

Perfect! We're only seeing real contig names that (usually) represent aligned reads. Let's count them by piping to wc -l (and omitting omit head of course – we want to count everything).

    grep -P '^M00538' brca_pairedend.sam | cut -f 3 | grep -v '*' | wc -l


#### Exercise: About how many records represent aligned sequences? What alignment rate does this represent?


### Use sort and uniq to create a histogram of mapping qualities

The mapping quality score is in field 5 of the brca_pairedend.sam file. We can do this to pull out only that field:

    grep -P '^M00538' brca_pairedend.sam | cut -f 5 | head
 

We will use the uniq create a histogram of these values. To create a histogram, we want to organize all equal mapping quality score lines into an adjacent block, then use uniq -c option to count them. The sort -n command does the sorting into blocks (-n means numerical sort). So putting it all together, and piping the output to a pager just in case, we get:

    grep -P '^M00538' brca_pairedend.sam | cut -f 5 | sort -n | uniq -c | more

We are looking at mapping quality values for both aligned and un-aligned records, but mapping quality only makes sense for aligned reads. This expression does not distinguish between mapping quality = 0 because the read mapped to multiple locations, and mapping quality = 0 because the sequence did not align. The proper solution will await the use of samtools to filter out unmapped reads.


After bowtie2 came out with a local alignment option, it wasn't long before bwa developed its own local alignment algorithm called **BWA-MEM (for Maximal Exact Matches)**, implemented by the bwa mem command. bwa mem has the following advantages:

- It incorporates a lot of the simplicity of using bwa with the complexities of local alignment, enabling straightforward alignment of datasets like the mirbase data we just examined
- It can align different portions of a read to different locations on the genome
- bwa mem was made available when we loaded the bwa module, so take a look at its usage information. The most important parameters, similar to those we've manipulated in the past two sections, are the following:

   -k	Controls the minimum seed length (default = 19)

   -w	Controls the "gap bandwidth", or the length of a maximum gap. This is particularly relevant for MEM, since it can determine whether a read is split into two separate alignments or is reported as one long alignment with a long gap in the middle (default = 100)

  -r	Controls how long an alignment must be relative to its seed before it is re-seeded to try to find a best-fit local match (default = 1.5, e.g. the value of -k multiplied by 1.5)

  -c	Controls how many matches a MEM must have in the genome before it is discarded (default = 10000)

  -t	Controls the number of threads to use


AS the documentation states it consists of three different algorithms: *BWA*, *BWA-SW* and *BWA-MEM*. The first algorithm, which is the oldest, is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA for 70-100bp Illumina reads.

##### Aligning with new BWA-MEM in both single-end (SE) and paired-end (PE) modes

BWA-MEM is the recommended algorithm to use now. You can check the options by executing:

    bwa mem

To align **PE** with BWA-MEM just execute the same command line with the two FASTQ files:

    bwa mem -t 2 -R "@RG\tID:foo\tSM:bar\tPL:Illumina\tPU:unit1\tLB:lib1"  data/ucsc.hg19 fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R1_001.fastq.gz  fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R2_001.fastq.gz > brca_pairedend_mem.sam 2> aln.brca_mem.log &
    
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

Now you can create the index by executing:

    bowtie2-build ucsc.hg19.fasta

Some files will be created in the ```ucsc.hg19``` folder, those files constitute the index that Bowtie2 uses.

**NOTE:** The index must created only once, it will be used for all the different alignments with Bowtie2.


##### Aligning in PE modes

Mapping in **PE** also requires only one execution:

    bowtie2 -q -p 4 -x ucsc.hg19.fasta -1 fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R1_001.fastq.gz -2  fq/510-7-ANGELA-SANTOS-BRCA_S8_L001_R2_001.fastq.gz -S brca_pairedend_bowtie.sam


### More exercises

- Work with your own data at home , [CQ_S16](https://dl.dropboxusercontent.com/u/1977573/sample1.zip).
- Try BWA-MEM algorithm vs Bowtie2 vs BWA-ALN
