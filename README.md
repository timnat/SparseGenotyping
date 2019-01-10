# SparseGenotyping
The SparseGenotyping program allows to resolve a consensus genotype for scaffolds based on low-coverage SNP calls by accumulating number of reads that support one or another genotype across the length of a scaffold.

This program was specifically written and used for the study presented in [1]. When sequencing depth at any position of reference genome is insufficient to reliably assess genotype for the position, instead, consensus genotypes can be defined for loci of some pre-determined or ambiguous length. In the current version program attempts to compute it for every scaffold of reference assembly.

In [1] the program was used to infer from low-coverage SNP calls whether each of individuals (from backcrossing F1 hybrid A. mexicanum/A.tigrinum to A. mexicanum) was homozygous A. mexicanum or heterozygous A. mexicanum/A. tigrinum. It is computed by summing the number of reads that could be assigned to A. mexicanum or A. tigrinum and assessing the resulting ratios (A. tigrinum reads / A. tigrinum + A. mexicanum reads) to define a consensus genotype for the scaffold. Individuals with ratios reflecting no or low representation of A. tigrinum alleles (<0.20) were scored as homozygous A. mexicanum, and individuals with ratios approximating 50% A. tigrinum alleles (0.25 - 0.75) were scored as heterozygous A. mexicanum/A. tigrinum. Individuals with intermediate frequencies or low coverage (fewer than 4 reads) were scored as missing genotypes. whether each individual was homozygous A. mexicanum or heterozygous A. mexicanum/A. tigrinum by summing the number of reads that could be assigned to A. mexicanum or A. tigrinum and assessing the resulting ratios (A. tigrinum reads / A. tigrinum + A. mexicanum reads) to define a consensus genotype for the scaffold. Individuals with ratios reflecting no or low representation of A. tigrinum alleles (<0.20) were scored as homozygous A. mexicanum, and individuals with ratios approximating 50% A. tigrinum alleles (0.25 - 0.75) were scored as heterozygous A. mexicanum/A. tigrinum. Individuals with intermediate frequencies or low coverage (fewer than 4 reads) were scored as missing genotypes.

Program is written on C++

INSTALLATION:

    download source code and compile with g++

    g++ SparseGenotyping.cpp -o SparseGenotyping

USAGE:
 
    ./SparseGenotyping [-r min_read_number [4]] [-m min_ratio [0.25]] [-M max_ratio [0.75]] [-e homoz_min_ratio [0.8]] vcf_file
    
    vcf_file - vcf file (output of GATK HaplotypeCaller) for NInd individuals [Default 50]. NInd is a predefined constant, that can be changed before compilation
    
    <min_read_number> - minimum number of reads to assign genotype [4]. At least min_read_number of reads, added accross all SNPs within a scaffold, is required to assign genotype
        
    <min_ratio> - minimum value of ratio to assign heterozygous genotype [0.25],  ratio = (number of reads of genotype1 / number of reads of genotype1 + number of reads of genotype2)
    
    <max_ratio> - maximum value of ratio to assign heterozygous genotype [0.75],  ratio = (number of reads of genotype1 / number of reads of genotype1 + number of reads of genotype2)
    
    <e> - if ratio > homoz_min_ratio, assign homozygous genotype1 [0.8]. Individuals with ratio < homoz_min_ratio and out of [min_ratio,max_ratio] interval intermediate frequencies are scored as missing genotypes

    Output:
    
    *.rd - file with read depth in a tabulated format with columns:
    
        Scaffold_name  
        
        Number_of_SNPsites  
        
        Number_of_reads_for_Indiv_with_homozigous_gt1
        
        Number_of_reads_for_Indiv_with_homozigous_gt2
        
        Number_of_reads_that_support_gt1_for_Indiv1
        
        Number_of_reads_that_support_gt2_for_Indiv1
        
        Number_of_reads_that_support_gt1_for_Indiv2
        
        Number_of_reads_that_support_gt2_for_Indiv2
        
        ...
        
    *.ratios - file that contains ratios=(Number_of_reads_that_support_gt1/Number_of_reads_that_support_gt2) for each individual
    
    *.gt - file that contains assigned genotypes
    
        0 - for homozygous genotype1
        
        1 - for heterozygous genotype, when reads from both alleles are observed almost same number of times
        
        ~ - can not assign genotype, either because of too few observed reads (even after summarizing across all SNPs) or when interpretation of ratio is rather uncertain, i.e. beyond [min_ratio;max_ratio].

EXAMPLE:

        Download test file AMex_Atigr_48indiv_test.vcf
        ./SparseGenotyping AMex_Atigr_48indiv_test.vcf -r 4 -m 0.3 -M 0.7 -e 0.9

CITATION:

1. A Chromosome-Scale Assembly of the Enormous (32 Gb) Axolotl Genome. Jeramiah J. Smith, Nataliya Timoshevskaya, Vladimir A. Timoshevskiy, Melissa C. Keinath, Drew Hardy, S. Randal Voss. www.biorxiv.org/content/biorxiv/early/2018/07/20/373548.full.pdf
