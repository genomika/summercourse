% [WorkShop de Verão Bioinformática](http://github.com/genomika/summercourse/)
% __Variant Calling and Post Processing__
% _(updated 09-02-2014)_


<!-- COMMON LINKS HERE -->

[SAMTools]: http://samtools.sourceforge.net/ "samtools"
[Picard]: http://picard.sourceforge.net/ "Picard"
[GATK]: http://www.broadinstitute.org/gatk/ "GATK"

Preliminaries
================================================================================

Software used in this practical:
--------------------------------

- [SAMTools] : SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [Picard] : Picard comprises Java-based command-line utilities that manipulate SAM files, and a Java API (SAM-JDK) for creating new programs that read and write SAM files.
- [GATK] : Genome Analysis Toolkit - A package to analyse next-generation re-sequencing data, primary focused on variant discovery and genotyping.


File formats explored:
----------------------

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf)
- [BAM](http://www.broadinstitute.org/igv/bam)
- VCF Variant Call Format: see [1000 Genomes](http://www.1000genomes.org/wiki/analysis/variant-call-format/vcf-variant-call-format-version-42) and [Wikipedia](http://en.wikipedia.org/wiki/Variant_Call_Format) specifications.


Exercise 1: Variant calling with paired-end data
================================================================================

Copy the necessary data in your working directory:

    mkdir -p $HOME/core_ngs/variant_calling/
    cp -r $HOME/core_ngs/align/*.bam $HOME/core_ngs/variant_calling/
    cp -r $HOME/core_ngs/align/*.bai $HOME/core_ngs/variant_calling/
    cd $HOME/core_ngs/variant_calling/

These datasets contain reads only for the *genes BRCA1, BRCA2**.


2. Prepare BAM file
--------------------------------------------------------------------------------

Go to the folder:

    cd $HOME/core_ngs/variant_calling/

The **read group** information is key for downstream GATK functionality. The GATK will not work without a read group tag. Make sure to enter as much metadata as you know about your data in the read group fields provided. For more information about all the possible fields in the @RG tag, take a look at the SAM specification.

    picard-tools AddOrReplaceReadGroups I=brca_pairedend_mem.bam O=fbrca_pairedend_mem.fixed.bam RGID=group1 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit1


We must sort and index the BAM file before processing it with Picard and GATK. To sort the bam file we use ``samtools``

    samtools sort brca_pairedend_mem.bam brca_pairedend_mem_sorted

Index the BAM file:

    samtools index brca_pairedend_mem_sorted.bam

**NOTE**: All those steps above are not recquired if you used **bwa mem** and **picard** as tools for alignment and sam to bam conversion respectively.


3. Mark duplicates (using Picard)
--------------------------------------------------------------------------------

Run the following **Picard** command to mark duplicates:

    picard-tools MarkDuplicates INPUT=brca_pairedend_mem.bam OUTPUT=brca_pairedend_mem.marked.bam METRICS_FILE=brca.txt CREATE_INDEX=true

This creates a sorted BAM file called ``brca_pairedend_mem.marked.bam`` with the same content as the input file, except that any duplicate reads are marked as such. It also produces a metrics file called ``brca.txt`` containing (can you guess?) metrics.

**QUESTION:** How many reads are removed as duplicates from the files (hint view the on-screen output from the two commands)?

4. Local realignment around INDELS (using GATK)
--------------------------------------------------------------------------------

There are 2 steps to the realignment process:

**First**, create a target list of intervals which need to be realigned
    
      java -Xmx2g -jar /usr/local/share/tools/GenomeAnalysisTK.jar -nt 2  -T RealignerTargetCreator  -R /mnt/data/ucsc.hg19.fasta -o brca_pairedend_mem.marked.bam.list -I brca_pairedend_mem.marked.bam -known  /mnt/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -known /mnt/data/1000G_phase1.indels.hg19.vcf -L brca.list


**Second**, perform realignment of the target intervals

    java -Xmx10g -jar /usr/local/share/tools/GenomeAnalysisTK.jar -I brca_pairedend_mem.marked.bam -R /mnt/data/ucsc.hg19.fasta -T IndelRealigner -targetIntervals brca_pairedend_mem.marked.bam.list -known /mnt/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -known /mnt/data/1000G_phase1.indels.hg19.vcf -o brca_pairedend_mem.marked.realigned.bam

This creates a file called ``brca_pairedend_mem.marked.realigned.bam`` containing all the original reads, but with better local alignments in the regions that were realigned.


5. Base quality score recalibration (using GATK)
--------------------------------------------------------------------------------

Two steps:

**First**, analyse patterns of covariation in the sequence dataset

    java -Xmx10g -jar /usr/local/share/tools/GenomeAnalysisTK.jar -nct 2 -l INFO -T BaseRecalibrator -I brca_pairedend_mem.marked.realigned.bam -R /mnt/data/ucsc.hg19.fasta -L brca.list -knownSites /mnt/data/dbsnp_138.hg19.vcf  -knownSites /mnt/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites /mnt/data/1000G_phase1.indels.hg19.vcf -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -cov ReadGroupCovariate -o brca_pairedend_mem.marked.realigned.recal_data.table

This creates a GATKReport file called ``brca_pairedend_mem.marked.realigned.recal_data.table`` containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data.

It is imperative that you provide the program with a set of **known sites**, otherwise it will refuse to run. The known sites are used to build the covariation model and estimate empirical base qualities. For details on what to do if there are no known sites available for your organism of study, please see the online GATK documentation.

**Second**, apply the recalibration to your sequence data

    java -Xmx10g -jar /usr/local/share/tools/GenomeAnalysisTK.jar -nct 2 -T PrintReads -I brca_pairedend_mem.marked.realigned.bam  -R /mnt/data/ucsc.hg19.fasta -BQSR brca_pairedend_mem.marked.realigned.recal_data.table -o brca_pairedend_mem.marked.realigned.recal.bam
    
This creates a file called ``brca_pairedend_mem.marked.realigned.recal.bam`` containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag ``–emit_original_quals`` to the ``PrintReads`` command, in which case the original qualities will also be written in the file, tagged OQ.


6. Variant calling (using GATK - **HaplotypeGenotyper**)
--------------------------------------------------------------------------------

SNPs and INDELS are called using separate instructions.

**SNP calling**

    java -jar ../gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../genome/f000_chr21_ref_genome_sequence.fa -I 004-dna_chr21_100_hq_pe_sorted_noDup_realigned_recalibrated.bam -glm SNP -o 005-dna_chr21_100_he_pe_snps.vcf

**INDEL calling**

    java -jar ../gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../genome/f000_chr21_ref_genome_sequence.fa -I 004-dna_chr21_100_hq_pe_sorted_noDup_realigned_recalibrated.bam -glm INDEL -o 005-dna_chr21_100_hq_pe_indel.vcf


7. Introduce filters in the VCF file
--------------------------------------------------------------------------------

Example: filter SNPs with low confidence calling (QD < 12.0) and flag them as "LowConf".

    java -jar ../gatk/GenomeAnalysisTK.jar -T VariantFiltration -R ../genome/f000_chr21_ref_genome_sequence.fa -V 005-dna_chr21_100_he_pe_snps.vcf --filterExpression "QD < 12.0" --filterName "LowConf" -o 006-dna_chr21_100_he_pe_snps_filtered.vcf

The command ``--filterExpression`` will read the INFO field and check whether variants satisfy the requirement. If a variant does not satisfy your filter expression, the field FILTER will be filled with the indicated ``--filterName``. These commands can be called several times indicating different filtering expression (i.e: --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2").

**QUESTION:** How many "LowConf" variants have we obtained?

    grep LowConf 006-dna_chr21_100_he_pe_snps_filtered.vcf | wc -l

