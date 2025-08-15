#rule all:
 #   input:
  #      expand("results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config["IDS_0"]),
   #     expand("results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config["IDS_1"]),
    #    expand("results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config["IDS_2"]),
     #   expand("results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config["IDS_3"]),
      #  expand("results/eQTL/sample_check/{sample}.mbv_output.txt", sample = config["IDS_4"])


rule sample_check_1:
    input:
        vcf = multiext("work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered", ".vcf.gz", ".vcf.gz.tbi"),
        bam = "work/RNA_seq/align/T0/{sample}.bam",
        bai = "work/RNA_seq/align/T0/{sample}.bam.bai"
    output:
        check = "results/eQTL/sample_check/{sample}.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''
rule sample_check_2:
    input:
        vcf = multiext("work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered", ".vcf.gz", ".vcf.gz.tbi"),
        bam = "work/RNA_seq/align/T1/{sample}.bam",
        bai = "work/RNA_seq/align/T1/{sample}.bam.bai"
    output:
        check = "results/eQTL/sample_check/{sample}.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''
rule sample_check_3:
    input:
        vcf = multiext("work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered", ".vcf.gz", ".vcf.gz.tbi"),
        bam = "work/RNA_seq/align/T2/{sample}.bam",
        bai = "work/RNA_seq/align/T2/{sample}.bam.bai"
    output:
        check = "results/eQTL/sample_check/{sample}.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''

rule sample_check_4:
    input:
        vcf = multiext("work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered", ".vcf.gz", ".vcf.gz.tbi"),
        bam = "work/RNA_seq/align/T3/{sample}.bam",
        bai = "work/RNA_seq/align/T3/{sample}.bam.bai"
    output:
        check = "results/eQTL/sample_check/{sample}.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''

rule sample_check_5:
    input:
        vcf = multiext("work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered", ".vcf.gz", ".vcf.gz.tbi"),
        bam = "work/RNA_seq/align/T4/{sample}.bam",
        bai = "work/RNA_seq/align/T4/{sample}.bam.bai"
    output:
        check = "results/eQTL/sample_check/{sample}.mbv_output.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam} --filter-mapping-quality 150 --filter-base-quality 20  --out {output.check} 
        '''