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
     

###Exercise: Generate "yeast_chip.bam", a sorted bam file named "yeast_chip_sort.bam", and an index file of the sorted bam file

    samtools view -bS yeast_chip.sam > yeast_chip.bam
    samtools sort yeast_chip.bam yeast_chip_sort
    samtools index yeast_chip_sort.bam
    
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
    
