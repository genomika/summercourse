[WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Alinhamento__
% _(updated 29-01-2014)_

<!-- COMMON LINKS HERE -->

[SAMTools]: http://samtools.sourceforge.net/ "SAMtools"
[Picard]: http://broadinstitute.github.io/picard/ "Picard"

# Output files of mapping

- Mappers spit out SAM/BAM files as output files. BAM is just a binary format of SAM.
- Here is a specification of SAM format SAM specification.
- Commonly, SAM files are processed in this order:
    - SAM files are converted into BAM files (samstools view)
    - BAM files are sorted by reference coordinates (samtools sort)
    - Sorted BAM files are indexed (samtools index)
    - Each step above can be done with commands below

Each step above can be done with commands below:

     samtools view -bS <samfile> > <bamfile>
     samtools sort <bamfile> <prefix of sorted bamfile>
     samtools index <sorted bamfile>
     

###Exercise: Generate "brca_pairedend.sort.bam", a sorted bam file named "brca_paired_sort.bam", and an index file of the sorted bam file

    cd /$HOME/core_ngs/align
    samtools view -bS brca_pairedend.sam > brca_pairedend.bam
    samtools sort brca_paired.bam brca_paired_sort
    samtools index brca_paired_sort.bam
    
### Basic mapping stats

    samtools flagstat <bamfile>

It prints out a host of useful stats about your mapping results, such as reads mapped. It can print a lot more information like % properly paired, # of duplicates, but it's simply relying on the information encoded in the second field of the SAM file - the bitwise flag field.

If you're using a mapper that doesn't generate a SAM file, stop what you are doing; go back and get a mapper that makes a SAM file.

####Exercise: Print out basic stats about yeast_chip.bam

    samtools idxstats <indexed bam file>

reports on stats related to the chromosome-based indexing done by samtools index. For each sequence of the reference, it provides:
- Sequence name (usually "chr1", etc.)
- BP in that sequence
- Reads mapping to that sequence
- Reads not mapping to that sequence    
    
### Additional SAMtools tricks

#### Extract/print sub alignments in BAM format.

If no region is specified in samtools view command, all the alignments will be printed; otherwise only alignments overlapping the specified regions will be output. A region can be presented, for example, in the following format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’

#### Count the number of mapped reads overlapping between chromossome 17 41197694 and chromosome 17 41197819

    samtools view brca_paired_sort.bam chr17:41197694-41197819 | wc -l

#### Extract/print mapped sub alignments in BAM format

A bam/sam file includes or does not include unmapped reads depending on mappers or options on mappers. If you use bwa with default options, the output bam includes unmapped reads. In order to extract mapped reads from a bam file, use -F option in samtools view command. -F INT means "Skip alignments with bits present in INT". In other words, -F INT filters reads that have the INT in their flag. Please take a look at page 4 on SAM specification. 0X4 flag is for segment unmapped.


#### Count the number of reads mapped on chromosome 17

    samtools view -F 0X04 brca_paired_sort.bam chr17 | wc -l
    
#### Extract/print reversely mapped sub alignments in BAM format

If you have a strand-specific RNA-seq data, a mapping direction of a read is critical. For example, you may want to count only reversely-mapped reads on a (-)-direction gene. Directionality of mapped read is also recorded on flag.

#### Count the number of reversely-mapped reads overlapping between chromosome 17: 41197694 and chromosome 17: 41197819
**Hint:** flag 0x10 = SEQ being reverse complemented

    samtools view -F 0X04 -f 0X10 brca_paired_sort.bam chr17:41197694-41197819 | wc -l

#### Count the total number of mappped reads in 'brca_pairedend.bam' 

    samtools flagstat brca_pairedend.bam
    samtools view -F 0x04 brca_paired_sort.bam | wc -l

# Output files of mapping (Alternative approach)

Besides samtools, there's another popular tool available written in java called [Picard](broadinstitute.github.io/picard/). Picard offers many options for manipulating or viewing SAM and BAM.

Commonly, SAM files are processed in this order:
    - SAM files are converted into BAM files (SortSam.jar)
    - BAM files are sorted by reference coordinates (SortSam.jar)
    - Sorted BAM files are indexed (SortSam.jar)
    - Each step above can be done with commands below

Each step above can be done with only one command below:

    java -Xmx4g -Djava.io.tmpdir=/tmp \
    -jar picard/SortSam.jar \
    SO=coordinate \
    INPUT=input.sam \
    OUTPUT=output.bam \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

For having a better overview I put the arguments in separate lines and join them in a script by the “\” sign.

The file types are recognized by their endings, so be sure to have them named as “.sam” and “.bam” files. The SO argument specifies the sort order which is by coordinate in our case. By setting the validation stringency to lenient, picard ignores some validation errors which frequently occur at the alignment step. 

By setting the CREATE_INDEX argument to true we automatically create an index file for the generated bam file. 

This file got the same name as the bam file but with the additional extension “.bai” (so sample1.bam got the index file sample1.bam.bai).

Let's do with our sample1

    picard-tools SortSam \
    SO=coordinate \
    INPUT=brca_pairedend.sam \
    OUTPUT=brca_pairedend_picard.bam \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

###Exercise: Now try the samtools and picard with the [CQ_S16](https://dl.dropboxusercontent.com/u/1977573/sample1.zip) sample!
