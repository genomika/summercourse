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
