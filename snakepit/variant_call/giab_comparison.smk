# Script to compare proportion of correct allele calls for GIAB samples and this dataset

ids = [1,5,6,7]

#rule all:
 #   input:
  #      expand('work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt', id = ids),
   #     'work/RNA_seq/samplecheck/variantcomparison/DV_nerged_raw_calls.txt',
    #    'work/RNA_seq/samplecheck/variantcomparison/DV_raw_GIAB_concordance.pdf'


rule GIAB:
    input:
        HG = 'data/GIAB/HG00{id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
       
    output:
        HG_calls = 'work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt'

    shell:
        '''
        software/bcftools-1.15.1/bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' {input.HG} > {output.HG_calls}
        '''


rule deepvariant_calls:
    input:
        DV = 'work/RNA_seq/variant_call/mergedcall/all.Unrevised.vcf.gz'

    output:
        DV_calls = 'work/RNA_seq/samplecheck/variantcomparison/DV_merged_raw_calls.txt'

    shell:
        '''
        software/bcftools-1.15.1/bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' {input.DV} > {output.DV_calls}
        '''
rule plotting_GIAB:
    input:
        DV_calls = 'work/RNA_seq/samplecheck/variantcomparison/DV_merged_raw_calls.txt',
        HG_calls = expand('work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt', id = ids),
        script = "scripts/GIAB_plotting.R"
    output:
        plot = 'work/RNA_seq/samplecheck/variantcomparison/DV_raw_GIAB_concordance.pdf',
        triallelic = "work/DNA_seq/imputation/variantcomparison/DV_no_triallelic_GIAB_concordance.pdf"
    shell:
        '''
        Rscript {input.script} {input.DV_calls} {input.HG_calls} {output.plot} {output.triallelic}
        '''
