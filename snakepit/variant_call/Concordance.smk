#rule all:
 #   input:
  #      expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{timepoint}_Genotypes.GT.FORMAT', timepoint=config["TIME"]),
      #  '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T1.pdf',
       # '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T2.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T3.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T4.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T0.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T2.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T3.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T4.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T0.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T1.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T3.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T4.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T0.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T1.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T2.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T4.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T0.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T1.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T2.pdf',
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T3.pdf'

# To get other files, just change the Ts to correspond to comparisons you are looking for and update rule all accordingly

rule GTExtract:
    input:
        vcfs = lambda wildcards: expand(f'/home/workspace/jogrady/heQTL/work/RNA_seq/variant_call/{config["TIME"][wildcards.timepoint]}/all.Unrevised.vcf.gz'),
        rename = expand('/home/workspace/jogrady/heQTL/data/{{timepoint}}/Rename_{{timepoint}}.txt'),
    output:
        GTs = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{{timepoint}}_Genotypes.GT.FORMAT')
    shell:
        '''
        /home/workspace/jogrady/heQTL/software/bcftools-1.15.1/bcftools reheader {input.vcfs} -s {input.rename} | bcftools view  -S /home/workspace/jogrady/heQTL/data/T0/sampleorder.txt -Ov | bcftools +fill-tags -Ov | bcftools view  -i  'MIN(FMT/DP)>5' -Ov |  bcftools view -i 'MAF > 0.1' -Ov | bcftools filter -e 'HWE < 0.000001' -Ov | bcftools filter -e 'F_MISSING > 0.05' -Ov | vcftools --vcf - --extract-FORMAT-info GT --out /home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{wildcards.timepoint}_Genotypes
        '''

rule ConcorT0:
    input:
        T0 = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/T0_Genotypes.GT.FORMAT',
        T1 = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{time}_Genotypes.GT.FORMAT', time = ["T1", "T2", "T3", "T4"]),
        script = '/home/workspace/jogrady/heQTL/scripts/Concordance_pdf.R'
    output:
        heatmap = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T0_V_T1.pdf',
    shell:
        '''
        Rscript {input.script} {input.T0} {input.T1}[0] T0 T1
        Rscript {input.script} {input.T0} {input.T1}[1] T0 T2
        Rscript {input.script} {input.T0} {input.T1}[2] T0 T3
        Rscript {input.script} {input.T0} {input.T1}[3] T0 T4
        '''
rule ConcorT1:
    input:
        T1 = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/T1_Genotypes.GT.FORMAT',
        T2 = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{time}_Genotypes.GT.FORMAT', time = ["T0", "T2", "T3", "T4"]),
        script = '/home/workspace/jogrady/heQTL/scripts/Concordance_pdf.R'
    output:
        heatmap = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T1_V_T0.pdf',
    shell:
        '''
        Rscript {input.script} {input.T1} {input.T2}[0] T1 T0
        Rscript {input.script} {input.T1} {input.T2}[1] T1 T2
        Rscript {input.script} {input.T1} {input.T2}[2] T1 T3
        Rscript {input.script} {input.T1} {input.T2}[3] T1 T4
        '''
rule ConcorT2:
    input:
        T2 = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/T2_Genotypes.GT.FORMAT',
        T3 = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{time}_Genotypes.GT.FORMAT', time = ["T0", "T1", "T3", "T4"]),
        script = '/home/workspace/jogrady/heQTL/scripts/Concordance_pdf.R'
    output:
        heatmap = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T2_V_T0.pdf'
    shell:
        '''
        Rscript {input.script} {input.T2} {input.T3}[0] T2 T0
        Rscript {input.script} {input.T2} {input.T3}[1] T2 T1
        Rscript {input.script} {input.T2} {input.T3}[2] T2 T3
        Rscript {input.script} {input.T2} {input.T3}[3] T2 T4
        '''
rule ConcorT3:
    input:
        T0 = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/T3_Genotypes.GT.FORMAT',
        T3 = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{time}_Genotypes.GT.FORMAT', time = ["T0", "T1", "T2", "T4"]),
        script = '/home/workspace/jogrady/heQTL/scripts/Concordance_pdf.R'
    output:
        heatmap = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T3_V_T0.pdf',
    shell:
        '''
        Rscript {input.script} {input.T0} {input.T3}[0] T3 T0
        Rscript {input.script} {input.T0} {input.T3}[0] T3 T1
        Rscript {input.script} {input.T0} {input.T3}[0] T3 T2
        Rscript {input.script} {input.T0} {input.T3}[0] T3 T4
        '''



rule ConcorT4:
    input:
        T0 = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/T4_Genotypes.GT.FORMAT',
        T5 = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/{time}_Genotypes.GT.FORMAT', time = ["T0", "T1", "T2", "T3"]),
        script = '/home/workspace/jogrady/heQTL/scripts/Concordance_pdf.R'
    output:
        heatmap = '/home/workspace/jogrady/heQTL/work/RNA_seq/samplecheck/Concordance_T4_V_T0.pdf',
    shell:
        '''
        Rscript {input.script} {input.T0} {input.T5}[1] T4 T0
        Rscript {input.script} {input.T0} {input.T5}[2] T4 T1
        Rscript {input.script} {input.T0} {input.T5}[3] T4 T2
        Rscript {input.script} {input.T0} {input.T5}[4] T4 T3
        '''

