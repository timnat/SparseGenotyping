# SparseGenotyping
SparseGenotyping - is a program that allows to define a consensus genotype for scaffolds based on low-coverage SNP calls by accumulating number of reads that support one or another genotype across the length of a scaffold.

This program was specifically written and used for the study presented in [https://www.biorxiv.org/content/biorxiv/early/2018/07/20/373548.full.pdf]. When sequencing depth at any position of reference genome is insufficient to reliably assess genotype for the position, instead, consensus genotypes can be defined for loci of some pre-determined or ambiguous length. In the current version program attempts to compute it for every scaffold of reference assembly.

In [] the program was used to infer from low-coverage SNP calls whether each of individuals (from backcrossing F1 hybrid A. mexicanum/A.tigrinum to A. mexicanum) was homozygous A. mexicanum or heterozygous A. mexicanum/A. tigrinum. It is computed by summing the number of reads that could be assigned to A. mexicanum or A. tigrinum and assessing the resulting ratios (A. tigrinum reads / A. tigrinum + A. mexicanum reads) to define a consensus genotype for the scaffold. Individuals with ratios reflecting no or low representation of A. tigrinum alleles (<0.20) were scored as homozygous A. mexicanum, and individuals with ratios approximating 50% A. tigrinum alleles (0.25 - 0.75) were scored as heterozygous A. mexicanum/A. tigrinum. Individuals with intermediate frequencies or low coverage (fewer than 4 reads) were scored as missing genotypes. whether each individual was homozygous A. mexicanum or heterozygous A. mexicanum/A. tigrinum by summing the number of reads that could be assigned to A. mexicanum or A. tigrinum and assessing the resulting ratios (A. tigrinum reads / A. tigrinum + A. mexicanum reads) to define a consensus genotype for the scaffold. Individuals with ratios reflecting no or low representation of A. tigrinum alleles (<0.20) were scored as homozygous A. mexicanum, and individuals with ratios approximating 50% A. tigrinum alleles (0.25 - 0.75) were scored as heterozygous A. mexicanum/A. tigrinum. Individuals with intermediate frequencies or low coverage (fewer than 4 reads) were scored as missing genotypes.

Program is written on C++

INSTALLATION:
    download source code and compile with g++
    g++ SparseGenotyping.cpp -o SparseGenotyping

USAGE:
    Input: should be in order
    SparseGenotyping <input.vcf> <min_read_number> <min_ratio> <max_ratio> <e>
    input.vcf - vcf file (output of GATK HaplotypeCaller), in which first and second individual (columns) are individuals with known genotypes, and genotype have to be determined for all successive individuals [Default 50]
    <min_read_number> - minimum number of reads to assign genotype [4]
    <min_ratio> - minimum value of ratio = (number of reads of genotype 1 / number of reads of gt1 + number of reads of genotype 2) to assign heterozygous genotype  [0.25]
    <max_ratio> - maximum value of ratio = (number of reads of genotype 1 / number of reads of gt1 + number of reads of genotype 2) to assign heterozygous genotype [0.75]
    <e> - if ratio > e, assign homozygous genotype 1 [0.8]. Individuals with intermediate frequencies or low coverage (fewer than mean_read_number reads) were scored as missing genotypes

    Output:
    <input.vcf>.rd - file with read depth in a tabulated format with columns:
        Scaffold_name  
        Number_of_SNPsites  
        Number_of_reads_for_Indiv_with_homozigous_gt1
        Number_of_reads_for_Indiv_with_homozigous_gt2
        Number_of_reads_that_support_gt1_for_Indiv1
        Number_of_reads_that_support_gt2_for_Indiv1
        Number_of_reads_that_support_gt1_for_Indiv2
        Number_of_reads_that_support_gt2_for_Indiv2
        ...
    <input.vcf>.ratios - file that contains ratios=(Number_of_reads_that_support_gt1/Number_of_reads_that_support_gt2) for each individual
    <input.vcf>.gt - file that contains assigned genotypes
        0 - for homozygous genotype 1
        1 - for heterozygous genotype, when both allele for gt1 and allele for gt2 observed almost same number of times
        ~ - can not assign genotype, either because of too few observed reads (even after summarizing across all SNPs) or when interpretation of ratio is rather uncertain, i.e. beyond [min_ratio;max_ratio].


