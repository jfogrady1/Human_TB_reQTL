# Script to compare proportion of correct allele calls for GIAB samples and this dataset

ids = [1,5,6,7]

#rule all:
 #   input:
  #      expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt', id = ids),
   #     '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_nerged_raw_calls.txt',
    #    '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_raw_GIAB_concordance.pdf'


rule GIAB:
    input:
        HG = '/home/workspace/jogrady/heQTL/data/GIAB/HG00{id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
       
    output:
        HG_calls = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt'

    shell:
        '''
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' {input.HG} > {output.HG_calls}
        '''


rule deepvariant_calls:
    input:
        DV = '/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/mergedcall/all.Unrevised.vcf.gz'

    output:
        DV_calls = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_merged_raw_calls.txt'

    shell:
        '''
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' {input.DV} > {output.DV_calls}
        '''
rule plotting_GIAB:
    input:
        DV_calls = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_merged_raw_calls.txt',
        HG_calls = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/HG00{id}_calls.txt', id = ids),
        script = "/home/workspace/jogrady/heQTL/scripts/GIAB_plotting.R"
    output:
        plot = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/variantcomparison/DV_raw_GIAB_concordance.pdf',
        triallelic = "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/DV_no_triallelic_GIAB_concordance.pdf"
    shell:
        '''
        Rscript {input.script} {input.DV_calls} {input.HG_calls} {output.plot} {output.triallelic}
        '''
