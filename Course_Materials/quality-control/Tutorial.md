% [WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Controle de Qualidade__
% _(updated 23-01-2014)_


<!-- Common URLs: Tools -->

[fastqc]:http://www.bioinformatics.babraham.ac.uk/projects/fastqc "FastQC home page"
[cutadapt]:http://code.google.com/p/cutadapt "cutadapt home page"

<!-- Common URLs: File Formats -->

[fastq-format-wikipedia]:http://en.wikipedia.org/wiki/FASTQ_format  "FastQC in Wikipedia"
[fastq-format-nar]:http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217 "FastQC original NAR publication"

<!-- External URLs -->

[mirbase]:http://www.mirbase.org/ "miRBase: a searchable database of published miRNA"


Preliminares
================================================================================

Softwares usados neste tutorial
----------------------------------------

- [FastQC] : A quality control tool for high-throughput sequence data.
- [cutadapt] : A tool to remove adapter sequences from high-throughput sequencing data.


### Setup
##### Special server login
Before we start, log into server like you did yesterday, but use this special hostname:

     172.16.225.177


File formats explored in this practical:
----------------------------------------

__FastQ__. See: 

- [Wikipedia][fastq-format-wikipedia].
- [NAR 2010][fastq-format-nar].


#### Data staging

Let's now set ourselves up to process this data in the server, using some of best practices for organizing our workflow.

    # Create a work area to work on data for this course,
    # with a sub-directory for pre-processing raw fastq files
    mkdir -p $HOME/core_ngs/fastq_prep
    # Make a symbolic links to the original yeast data in $WORK.
    cd $HOME/core_ngs/fastq_prep
    ln -s -f $HOME/work/Sample_Yeast_L005_R1.cat.fastq.gz
    ln -s -f $HOME/work/Sample_Yeast_L005_R2.cat.fastq.gz


Overview
================================================================================

1. Use [FastQC] to explore the raw data.
1. Use [cutadapt] to remove adapters.
1. Use [cutadapt] to filter reads based on quality.
1. Use [FastQC] to explore the filtered data.


Explore the raw data using some Linux shell commands
--------------------------------------------------------------------------------

##### Illumina sequence data format (FASTQ)

Paired-end sequencing gives you paired end sequencing data in two matching fastq format files, containing reads for each end sequenced. See where your data really is and how big it is.

    # the -l options says "long listing" which shows where the link goes,
    # but doesn't show details of the real file
    ls -l
    # the -L option says to follow the link to the real file, -l = long listing (includes size)
    # and -h says "human readable" (e.g. MB, GB)
    ls -lLh

#### 4-line FASTQ format

Each read end sequenced is represented by a 4-line entry in the fastq file that looks like this:

    @HWI-ST1097:127:C0W5VACXX:5:1101:4820:2124 1:N:0:CTCAGA
    TCTCTAGTTTCGATAGATTGCTGATTTGTTTCCTGGTCATTAGTCTCCGTATTTATATTATTATCCTGAGCATCATTGATGGCTGCAGGAGGAGCATTCTC
    +
    CCCFFFFDHHHHHGGGHIJIJJIHHJJJHHIGHHIJFHGICA91CGIGB?9EF9DDBFCGIGGIIIID>DCHGCEDH@C@DC?3AB?@B;AB??;=A>3;; 

**Line 1** is the unique read name. The format is as follows, using the read name above:

     machine_id:lane:flowcell_grid_coordinates  end_number:failed_qc:0:barcode
     @HWI-ST1097:127:C0W5VACXX:5:1101:4820:2124 1:N:0:CTCAGA

The line as a whole will be unique for this read fragment. However, the corresponding R1 and R2 reads will have identical machine_id:lane:flowcell_grid_coordinates information. This common part of the name ties the two read ends together.

Most sequencing facilities will not give you qc-failed reads (failed_qc = Y) unless you ask for them.

**Line 2** is the sequence reported by the machine, starting with the first base of the insert (the 5' adapter has been removed by the sequencing facility). These are ACGT or N uppercase characters.

**Line 3** is always '+' from GSAF (it can optionally include a sequence description)

**Line 4** is a string of Ascii-encoded base quality scores, one character per base in the sequence.
For each base, an integer Phred-type quality score is calculated as integer score = -10 log(probabilty base is wrong) then added to 33 to make a number in the Ascii printable character range. As you can see from the table below, alphabetical letters - good, numbers – ok, most special characters – bad (except :;<=>?@).

![Terminal Example](https://wikis.utexas.edu/download/attachments/66696890/ascii_qualities.png?version=1&modificationDate=1400197837000&api=v2)


See the Wikipedia FASTQ format page for more information.

#### Exercise: What character in the quality score string in the fastq entry above represents the best base quality?


#### gzip and gunzip

With no options, gzip compresses the file you give it in-place. Once all the content has been compressed, the original uncompressed file is removed, leaving only the compressed version (the original file name plus a .gz extension). The gzunzip function works in a similar manner, except that its input is a compressed file with a .gz file and produces an uncompressed file without the .gz extension.

    # make sure you're in your $HOME/core_ngs/fastq_prep directory
    cp $HOME/data/small.fq .
    # check the size, then compress it in-place
    ls -lh
    gzip small.fq
    # check the compressed file size, then uncompress it
    ls -lh 
    gunzip small.fq.gz

**NOTE** Both gzip and gunzip are extremely I/O intensive when run on large files. While bioinformatic servers has tremendous compute resources and  file system is not so great, it has its limitations. The intensity of compression/decompression operations is another reason you should compress your sequencing files once (if they aren't already) then leave them that way.

### About compressed files

Sequencing data files can be very large - from a few megabytes to gigabytes. And with NGS giving us longer reads and deeper sequencing at decreasing price points, it's not hard to run out of storage space. As a result, most sequencing facilities will give you compressed sequencing data files. The most common compression program used for individual files is gzip and its counterpart gunzip whose compressed files have the .gz extension. The tar and zip programs are most commonly used for compressing directories.

Let's take a look at the size difference between uncompressed and compressed files. We use the -l option of ls to get a long listing that includes the file size, and -h to have that size displayed in "human readable" form rather than in raw byte sizes.

    ls -lh $HOME/data/yeast_stuff/*.fastq
    ls -lh $HOME/data/yeast_stuff/*.fastq.gz

#### Exercise: About how big are the compressed files? The uncompressed files? About what is the compression factor?

You may be tempted to want to uncompress your sequencing files in order to manipulate them more directly – but resist that temptation. Nearly all modern bioinformatics tools are able to work on .gz files, and there are tools and techniques for working with the contents of compressed files without ever uncompressing them.

#### head and tail, more or less

One of the challenges of dealing with large data files, whether compressed or not, is finding your way around the data – finding and looking at relevant pieces of it. Except for the smallest of files, you can't open them up in a text editor because those programs read the whole file into memory, so will choke on sequencing data files! Instead we use various techniques to look at pieces of the files at a time.

The first technique is the use of "pagers" – we've already seen this with the **more** command. Review its use now on our small uncompressed file:

    # Use spacebar to advance a page; Ctrl-c to exit
    more small.fq

Another pager, with additional features, is **less**. The most useful feature of less is the ability to search – but it still doesn't load the whole file into memory, so searching a really big file can be slow.

Here's a summary of the most common less navigation commands, once the less pager is active. It has tons of other options (try less --help).
- q – quit
- Ctrl-f or space – page forward
- Ctrl-b – page backward
- /<pattern> – search for <pattern> in forward direction
- n – next match
- N – previous match
- ?<pattern> – search for <pattern> in backward direction
- n – previous match going back
- N – next match going forward
- If you start less with the -N option, it will display line numbers.

### Exercise: What line of small.fq contains the read name with grid coordinates 2316:10009:100563?

For a really quick peak at the first few lines of your data, there's nothing like the **head** command. By default it displays the first 10 lines of data from the file you give it or from its standard input. With an argument -NNN (that is a dash followed by some number), it will show that many lines of data.

    # shows 1st 10 lines
    head small.fq
    # shows 1st 100 lines -- might want to pipe this to more to see a bit at a time
    head -100 small.fq | more

The yang to head's ying is **tail**, which by default it displays the last 10 lines of its data, and also uses the -NNN syntax to show the last NNN lines. (Note that with very large files it may take a while for tail to start producing output because it has to read through the file sequentially to get to the end.)

See **piping** section for a diagram of what the piping operator ( | ) is doing.

But what's really cool about tail is its -n +NNN syntax. This displays all the lines starting at line NNN. Note this syntax: the -n option switch follows by a plus sign ( + ) in front of a number – the plus sign is what says "starting at this line"! Try these examples:

    # shows the last 10 lines
    tail small.fq
    # shows the last 100 lines -- might want to pipe this to more to see a bit at a time
    tail -100 small.fq | more
    # shows all the lines starting at line 900 -- better pipe it to a pager!
    tail -n +900 small.fq | more
    # shows 15 lines starting at line 900 because we pipe to head -15
    tail -n +900 small.fq | head -15

#### gunzip -c tricks

Ok, now you know how to navigate an uncompressed file using head and tail, more or less. But what if your FASTQ file has been compressed by gzip? You don't want to uncompress the file, remember? So you use the gunzip -c trick. This uncompresses the file, but instead of writing the uncompressed data to another file (without the .gz extension) it write it to its standard output where it can be piped to programs like your friends head and tail, more or less.

Let's illustrate this using one of the compressed files in your fq subdirectory:

    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | more
    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | head
    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | tail
    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | tail -n +900 | head -15
    # Note that less will display .gz file contents automatically
    less Sample_Yeast_L005_R1.cat.fastq.gz

**NOTE** There will be times when you forget to pipe your gunzip -c output somewhere – even the experienced among us still make this mistake! This leads to pages and pages of data spewing across your terminal. If you're lucky you can kill the output with Ctrl-c. But if that doesn't work (and often it doesn't) just close your Terminal window. This terminates the process on the server (like hanging up the phone), then you can log back in.

### Counting your sequences

One of the first thing to check is that your FASTQ files are the same length, and that length is evenly divisible by 4. The wc command (word count) using the -l switch to tell it to count lines, not words, is perfect for this. It's so handy that you'll end up using wc -l a lot to count things. It's especially powerful when used with filename wildcarding.

    wc -l small.fq
    head -100 small.fq > small2.fq
    wc -l small*.fq

You can also pipe the output of gunzip -c to wc -l to count lines in your compressed FASTQ file.

### Exercise: How many lines are in the Sample_Yeast_L005_R1.cat.fastq.gz file? How many sequences is this?

#### How to do math on the command line

The bash shell has a really strange syntax for arithmetic: it uses a double-parenthesis operator after the $ sign (which means evaluate this expression). Go figure.

    echo $((2368720 / 4))

Remember how the backticks evaluate the enclosed expression and substitute its output into a string? Here's how you would combine this math expression with gunzip -c line counting on your file using the magic of backtick evaluation:

    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | echo "$((`wc -l` / 4))"


### Processing multiple compressed files

You've probably figured out by now that you can't easily use filename wildcarding along with gunzip -c and piping to process multiple files. For this, you need to code a for loop in bash. Fortunately, this is pretty easy. Try this:

    for fname in *.gz; do
       echo "$fname has $((`gunzip -c $fname | wc -l` / 4)) sequences"
    done

Here fname is the name I gave the variable that is assigned a different file generated by the filename wildcarding list, each time through the loop. The actual file is then referenced as $file inside the loop.

Note the general structure of the for loop. Different portions of the structure can be separated on different lines (like <something> and <something else> below) or put on one line separated with a semicolon ( ; ) like before the do keyword below.

    for <variable name> in <expression>; do
       <something>
     <something else>
    done

You can also use high-level programming languages which can communicate with your os system. For instance, Python, an open-source language largely used by scientific community including bioinformatics. It is already native in several systems like debian, fedora, mac os, etc.  You can call a python script like the command below:

    python yourscript.py
    
This command will execute the code inside the yourscript.py .   So if we wanted to write a code that performs the same action as before, in this case, counting the number of sequences in each fastq.gz in a directory, the piece of code is presented as follows:

    import glob
    from subprocess import check_output


    all_fastqs = glob.glob('$HOME/core_ngs/fastq_prep/*.gz')
    for fastq in all_fastqs:
        print '%s has %d sequences' % (fastq, check_output(["gunzip -c %s | wc -l / 4" % fastq]) ) 
	
You may experiment by yourserf now:

    cd $HOME/core_ngs/fastq_prep
    ln -s -f $HOME/tools/preparer.py
    python preparer.py



Explore the raw data quality using FastQC
--------------------------------------------------------------------------------

The first order of business after receiving sequencing data should be to check your data quality. This often-overlooked step helps guide the manner in which you process the data, and can prevent many headaches.

**FastQC** is a tool that produces a quality analysis report on FASTQ files.

Useful links:
- FastQC [report](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) for a Good Illumina dataset
- FastQC [report](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) for a Bad Illumina dataset
- Online [documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) for each FastQC report

First and foremost, the FastQC "Summary" should generally be ignored. Its "grading scale" (green - good, yellow - warning, red - failed) incorporates assumptions for a particular kind of experiment, and is not applicable to most real-world data. Instead, look through the individual reports and evaluate them according to your experiment type.

The FastQC reports I find most useful, and why:
- Should I trim low quality bases?
	-  consult the Per base sequence quality report
- based on all sequences
	- Do I need to remove adapter sequences?
- consult the Overrepresented Sequences report
	- based on the 1st 200,000 sequences
- How complex is my library?
	- consult the Sequence Duplication Levels report

but remember that different experiment types are expected to have vastly different duplication profiles

**NOTE** For many of its reports, FastQC analyzes only the first 200,000 sequences in order to keep processing and memory requirements down. Consult the Online documentation for each FastQC report for full details.

**NOTE** Some of FastQC's graphs have a 1-100 vertical scale that is tricky to interpret. The 100 is a relative marker for the rest of the graph. For example, sequence duplication levels are relative to the number of unique sequences,

#### Running FastQC

It has a number of options (see fastqc --help | more) but can be run very simply with just a FASTQ file as its argument.

First create a directory to store the results of the fastqc analysis:

    mkdir f020_res_fastqc

Then execute `fastqc` storing the results in the created directory (option `-o`):

    fastqc -o f020_res_fastqc f010_raw_mirna.fastq

Find the results in the __fastqc_report.html__ file and discus them.


There are many _Overrepresented sequences_. 
Explore whether some of them correspond to miRNAs using the [miRBase search](http://www.mirbase.org/search.shtml) __By sequence__ utility.

### Exercise: What did FastQC create?

#### Looking at FastQC output
You can't run a web browser directly from your "dumb terminal" command line environment. The FastQC results have to be placed where a web browser can access them. We put a copy at this URL:

http://web.corral.tacc.utexas.edu/BioITeam/yeast_stuff/Sample_Yeast_L005_R1.cat_fastqc/fastqc_report.html

### Exercise: Based on this FastQC output, should we trim this data?


Handling adapters
--------------------------------------------------------------------------------

####Trimming sequences

There are two main reasons you may want to trim your sequences:
- As a quick way to remove 3' adapter contamination, when extra bases provide little additional information 
	- For example, 75+ bp ChIP-seq reads – 50 bases are more than enough for a good mapping, and trimming to 50 is easier than adapter removal, especially for paired end data.
	- You would not choose this approach for RNA-seq data, where 3' bases may map to a different exon, and that is valuable information. Instead you would specifically remove adapter sequences.
- Low quality base reads from the sequencer can cause an otherwise mappable sequence not to align
- This is more of an issue with sequencing for genome assembly – both bwa and bowtie2 seem to do fine with a few low quality bases, soft clipping them if necessary.

There are a number of open source tools that can trim off 3' bases and produce a FASTQ file of the trimmed reads to use as input to the alignment program.

#### FASTX Toolkit

The FASTX Toolkit provides a set of command line tools for manipulating both FASTA and FASTQ files. The available modules are described on their website. They include a fast fastx_trimmer utility for trimming FASTQ sequences (and quality score strings) before alignment.

Here's an example of how to run fastx_trimmer to trim all input sequences down to 50 bases. By default the program reads its input data from standard input and writes trimmed sequences to standard output:

    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | fastx_trimmer -l 50 -Q 33 > trim50_R1.fq

- The -l 50 option says that base 50 should be the last base (i.e., trim down to 50 bases)

- The -Q 33 option specifies how base qualities on the 4th line of each FASTQ entry are encoded. The FASTX Toolkit is an older program written in the time when Illumina base qualities were encoded differently, so its default does not work for modern FASTQ files. These days Illumina base qualities follow the Sanger FASTQ standard (Phred score + 33 to make an ASCII character).

#### Exercise: How would you tell fastx_trimmer to compress (gzip) its output file?


#### Adapter trimming with cutadapt

Data from RNA-seq or other library prep methods that result in short fragments can cause problems with moderately long (50-100bp) reads, since the 3' end of sequence can be read through to the 3' adapter at a variable position. This 3' adapter contamination can cause the "real" insert sequence not to align because the adapter sequence does not correspond to the bases at the 3' end of the reference genome sequence.

Unlike general fixed-length trimming (e.g. trimming 100 bp sequences to 50 bp), adapter trimming removes differing numbers of 3' bases depending on where the adapter sequence is found.

You must tell an adapter trimming program what your R1 and R2 adapters look like. The GSAF website describes the flavaors of Illumina adapter and barcode sequence in more detail https://wikis.utexas.edu/display/GSAF/Illumina+-+all+flavors.

The cutadapt program is an excellent tool for removing adapter contamination. The program is not available through TACC's module system but we linked a copy into your $HOME/local/bin directory.

A common application of cutadapt is to remove adapter contamination from small RNA library sequence data, so that's what we'll show here, assuming GSAF-format small RNA library adapters.

When you run cutadapt you give it the adapter sequence to trim, and this is different for R1 and R2 reads. Here's what the options look like (without running it on our files yet).

    cutadapt -m 22 -O 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    cutadapt -m 22 -O 4 -a TGATCGTCGGACTGTAGAACTCTGAACGTGTAGA

 The -m 22 option says to discard any sequence that is smaller than 22 bases (minimum) after trimming. This avoids problems trying to map very short, highly ambiguous sequences.

the -O 4 (Overlap) option says not to trim 3' adapter sequences unless at least the first  bases of the adapter are seen at the 3' end of the read. This prevents trimming short 3' sequences that just happen by chance to match the first few adapter sequence bases.

The -a argument to cutadapt is documented as the "sequence of adapter that was ligated to the 3' end". So we care about the <Read 2 primer> for R1 reads, and the <Read 1 primer> for R2 reads.

     cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC

There are 2 known adapters used in this experiment: 

    CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA
    CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC

Use the command grep to see whether they are still present in your data:

    grep "CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA" f010_raw_mirna.fastq 
	grep "CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC" f010_raw_mirna.fastq 

Do the sequences appear systematically at the beginning or at the end of the reads?



But the adapters could also appear in the _reverse_, _complementary_ or _reverse complementary_ mode.

Compute the _reverse_, _complementary_ and the _reverse complementary_ sequences of the two adapters,
and find out which of them appear in your data.

To compute those sequences you can use some online resources as the one in:  
<http://www.bioinformatics.org/sms/rev_comp.html>


<!--
Or to use R-Bioconductor to compute their reverse, complementary and reverse complementary.

library (Biostrings)
myseq <- DNAString ("CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC")
reverse (myseq)
complement (myseq)
reverseComplement (myseq)

-->

Use grep form Linux shell to find out which of the versions of the adapter is in your data.


### Adapter 1

    grep CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA f010_raw_mirna.fastq | wc -l  ## present in the sample (at the beginning of the reads)
	
    grep GACCCTTTAGTGGTATTTGCACTTTACAGAAACCTAAACCCTTAGAATATTCAAGACATACTCTGGTGAGATTTTT f010_raw_mirna.fastq | wc -l 
	
    grep TTTTTAGAGTGGTCTCATACAGAACTTATAAGATTCCCAAATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG f010_raw_mirna.fastq | wc -l  ## present in the sample (at the end of the read) ... but not so numerous
	
    grep AAAAATCTCACCAGAGTATGTCTTGAATATTCTAAGGGTTTAGGTTTCTGTAAAGTGCAAATACCACTAAAGGGTC f010_raw_mirna.fastq | wc -l 

But sometimes the adapter does not appear complete. 
It may be there just the first part:

    grep CTGGGAAATCACCATAAACGTGAAATGTCTTTGGA f010_raw_mirna.fastq | wc -l 
	
    grep GACCCTTTAGTGGTATTTGCACTTTACAGAAACCT f010_raw_mirna.fastq | wc -l 
	
    grep TTTTTAGAGTGGTCTCATACAGAACTTATAAGATT f010_raw_mirna.fastq | wc -l 
	
    grep AAAAATCTCACCAGAGTATGTCTTGAATATTCTAA f010_raw_mirna.fastq | wc -l 

or the end part of it:

    grep AATCTTATAAGTTCTGTATGAGACCACTCTAAAAA f010_raw_mirna.fastq | wc -l 
	
    grep TTAGAATATTCAAGACATACTCTGGTGAGATTTTT f010_raw_mirna.fastq | wc -l 
	
    grep TCCAAAGACATTTCACGTTTATGGTGATTTCCCAG f010_raw_mirna.fastq | wc -l 
	
    grep AGGTTTCTGTAAAGTGCAAATACCACTAAAGGGTC f010_raw_mirna.fastq | wc -l 

NOTE: in the code above I did cut just the 35 first or last nucleotides of the primer in its different arrangements, 
but this is an arbitrary length. 
We are just trying to discover which of the arrangements are present in our sample
and whether there are allocated in the 5' or in the 6' end.


### Adapter 2

    grep CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC f010_raw_mirna.fastq | wc -l   ## present in the sample (at the end of the read) ... but not so numerous
	
    grep GAAAAAAAGCAGGAAAGGTGTTCTATATATTTCGGTTCTTTAGCTTTATGAAAGTTCAATGCCATTCG f010_raw_mirna.fastq | wc -l 
	
    grep GCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG f010_raw_mirna.fastq | wc -l   ## present in the sample (at the beginning of the reads)
	
    grep CGAATGGCATTGAACTTTCATAAAGCTAAAGAACCGAAATATATAGAACACCTTTCCTGCTTTTTTTC f010_raw_mirna.fastq | wc -l 

As before, sometimes the adapter does not appear complete. 
It may be there just the first part:

    grep CTTTTTTTCGTCCTTTCCACAAGATATATA f010_raw_mirna.fastq | wc -l 
	
    grep GAAAAAAAGCAGGAAAGGTGTTCTATATAT f010_raw_mirna.fastq | wc -l 
	
    grep GCTTACCGTAACTTGAAAGTATTTCGATTT f010_raw_mirna.fastq | wc -l 
	
    grep CGAATGGCATTGAACTTTCATAAAGCTAAA f010_raw_mirna.fastq | wc -l 

or the end part of it:

    grep AAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC f010_raw_mirna.fastq | wc -l 
	
    grep TTCGGTTCTTTAGCTTTATGAAAGTTCAATGCCATTCG f010_raw_mirna.fastq | wc -l 
	
    grep CTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG f010_raw_mirna.fastq | wc -l 
	
    grep GAACCGAAATATATAGAACACCTTTCCTGCTTTTTTTC f010_raw_mirna.fastq | wc -l 



Use cutadapt to make an adapter trimming of the reads.
--------------------------------------------------------------------------------

Check the options:

- __-a__ for adapter to the 3' end.
- __-g__ for adapter to the 5' end.

you can find the help of the program typing `cutadapt -h` in the shell. 

To get read of the the adapters found in our data we run [cutadapt] several times:

    cutadapt -g CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA -o f030_mirna_trim1.fastq f010_raw_mirna.fastq
	
    cutadapt -a TTTTTAGAGTGGTCTCATACAGAACTTATAAGATTCCCAAATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG -o f030_mirna_trim2.fastq f030_mirna_trim1.fastq
	
    cutadapt -g GCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG         -o f030_mirna_trim3.fastq f030_mirna_trim2.fastq
	
    cutadapt -a CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC         -o f030_mirna_trim4.fastq f030_mirna_trim3.fastq


Now you can `grep` again searching for the adapters

### Adapter 1

    grep CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA f030_mirna_trim4.fastq | wc -l
	
    #grep GACCCTTTAGTGGTATTTGCACTTTACAGAAACCTAAACCCTTAGAATATTCAAGACATACTCTGGTGAGATTTTT f030_mirna_trim4.fastq | wc -l 
	
    grep TTTTTAGAGTGGTCTCATACAGAACTTATAAGATTCCCAAATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG f030_mirna_trim4.fastq | wc -l
	
    #grep AAAAATCTCACCAGAGTATGTCTTGAATATTCTAAGGGTTTAGGTTTCTGTAAAGTGCAAATACCACTAAAGGGTC f030_mirna_trim4.fastq | wc -l 

### Adapter 2

    grep CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC f030_mirna_trim4.fastq | wc -l
	
    #grep GAAAAAAAGCAGGAAAGGTGTTCTATATATTTCGGTTCTTTAGCTTTATGAAAGTTCAATGCCATTCG f030_mirna_trim4.fastq | wc -l 
	
    grep GCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG f030_mirna_trim4.fastq | wc -l
	
    #grep CGAATGGCATTGAACTTTCATAAAGCTAAAGAACCGAAATATATAGAACACCTTTCCTGCTTTTTTTC f030_mirna_trim4.fastq | wc -l 


The log file will look something like this:

    Command line parameters: -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC Sample_H54_miRNA_L004_R1.cat.fastq.gz
    Maximum error rate: 10.00%
       No. of adapters: 1
    Processed reads:      2001337
    Processed bases:    202135037 bp (202.1 Mbp)
    Trimmed reads:      1991662 (99.5%)
    Trimmed bases:    156404705 bp (156.4 Mbp) (77.38% of total)
    Too short reads:       246782 (12.3% of processed reads)
     Too long reads:            0 (0.0% of processed reads)
        Total time:     55.87 s
     Time per read:      0.028 ms

**Notes:**
- The Processed reads line tells you how many sequences were in the original FASTQ file.
- The Too short reads line tells you how may sequences were filtered out because they were shorter than our minimum length (20) after adapter removal (these may have ben primer dimers). 
- Here ~12% of the original sequences were removed, which is reasonable.
- Trimmed reads tells you how many of the reads you gave it had at least part of an adapter sequence that was trimmed.
- Here adapter was found in nearly all (99.5%) of the reads. This makes sense given this is a short (15-25 bp) RNA library.
- Trimmed bases tells you the total number of adapter bases that were removed, and its proportion of the total number of bases.
= Here > 3/4 of bases were trimmed. Again, this makes sense since the reads were 100 bm long.



Explore the quality of the trimmed file using FastQC
--------------------------------------------------------------------------------

Check the data quality again using fastqc:

	mkdir f040_res_fastqc_trimmed
	fastqc -o f040_res_fastqc_trimmed f030_mirna_trim4.fastq

Some of the reads seems to be too short and some others may not have enough quality. 


Use cutadapt to filter reads by quality and length. 
--------------------------------------------------------------------------------

Check the options:

- __-q__ quality cutoff.
- __-m__ minimum length.
- __-M__ maximum length.

you can find the help of the program typing `cutadapt -h` in the shell. 

Run cutadapt for length and quality purge of the reads.

    cutadapt -m 17 -M 30 -q 10 -o f040_mirna_cut.fastq f030_mirna_trim4.fastq


Check the data quality again using fastqc:

	mkdir f050_res_fastqc_trimmed_purged
	
	fastqc -o f050_res_fastqc_trimmed_purged f040_mirna_cut.fastq


Explore again the _Overrepresented sequences_ in [mirbase] (Search -> By sequence).

Count how many reads are left for the analysis (divide by 4)

    wc -l f010_raw_mirna.fastq
	
    wc -l f040_mirna_cut.fastq

### FASTQ statistics with samstat

The samstat program can also produce a quality report for FASTQ files, and it can also report on aligned sequences in a BAM file.

Again, this program is not available through the TACC module system but is available in Anna's work directory and has been linked into your $HOME/local/bin. You should be able just to type samstat and see some documentation.

Here's how you run samstat on a compressed FASTQ files. Let's only run it on a few sequences to avoid overloading the system:

    gunzip -c Sample_Yeast_L005_R1.cat.fastq.gz | head -100000 | samstat -f fastq -n samstat.fastq

This would produce a file named samstat.fastq.html which you need to view in a web browser. Here's some example output from bacterial sequences:

    http://web.corral.tacc.utexas.edu/BioITeam/SRR030257_1.fastq.html

