#rule all:
    #input:
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T0/{sample0}.bam', sample0 = config["fastqfiles_0"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/{sample1}.bam', sample1 = config["fastqfiles_1"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T2/{sample2}.bam', sample2 = config["fastqfiles_2"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T3/{sample3}.bam', sample3 = config["fastqfiles_3"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T4/{sample4}.bam', sample4 = config["fastqfiles_4"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_{time}.txt', time = config["TIME"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt', timepoint = config["TIME"]),
        #'/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_ALL.txt'

rule featureCounts0: 
    input:
        bam = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T0/{sample0}.bam', sample0 = config["fastqfiles_0"]),
        annotation="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_T0.txt'
    threads: 12
    shell:
        '''
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p --countReadPairs -C -T {threads} -s 2 -t gene -g gene_id
        '''

rule featureCounts1:
    input:
        bam = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T1/{sample1}.bam', sample1 = config["fastqfiles_1"]),
        annotation="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_T1.txt'
    threads: 12
    shell:
        '''
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p --countReadPairs -C -T {threads} -s 2 -t gene -g gene_id
        '''
rule featureCounts2:
    input:
        bam = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T2/{sample2}.bam', sample2 = config["fastqfiles_2"]),
        annotation="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_T2.txt'
    threads: 12
    shell:
        '''
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p --countReadPairs -C -T {threads} -s 2 -t gene -g gene_id
        '''
rule featureCounts3:
    input:
        bam = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T3/{sample3}.bam', sample3 = config["fastqfiles_3"]),
        annotation="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_T3.txt'
    threads: 12
    shell:
        '''
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p --countReadPairs -C -T {threads} -s 2 -t gene -g gene_id
        '''
rule featureCounts4:
    input:
        bam = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/T4/{sample4}.bam', sample4 = config["fastqfiles_4"]),
        annotation="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_T4.txt'
    threads: 12
    shell:
        '''
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p --countReadPairs -C -T {threads} -s 2 -t gene -g gene_id
        '''

rule cleanup:
    input:
        raw = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_{timepoint}.txt'

    output:
        cleaned = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_{timepoint}.txt'

    shell:
        r'''
        cut -f 1,7-55 {input.raw} > /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/cut_{wildcards.timepoint}.txt

        tail -n +2 /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/cut_{wildcards.timepoint}.txt > /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/tail_{wildcards.timepoint}.txt

        awk -F'\t' 'BEGIN {{OFS=FS}} {{sub("/home/workspace/jogrady/heQTL/work/RNA_seq/align/{wildcards.timepoint}/", "", $1); for (i=2; i<=NF; i++) sub("/home/workspace/jogrady/heQTL/work/RNA_seq/{wildcards.timepoint}/align/", "", $i)}} NF>0' /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/tail_{wildcards.timepoint}.txt > /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_{wildcards.timepoint}.txt
        rm /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/cut_{wildcards.timepoint}.txt

        Rscript /home/workspace/jogrady/heQTL/scripts/order_counts.R /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_{wildcards.timepoint}.txt {wildcards.timepoint}

        '''
rule featureCountsALL:
    input:
        bam = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{sample4}.bam', sample4 = config["merged_patient_IDS"]),
        annotation="/home/workspace/jogrady/heQTL/data/ref_genome/gencode.v43.annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_ALL.txt'
    threads: 20
    shell:
        '''
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p --countReadPairs -C -T 20 -s 2 -t gene -g gene_id
        '''



rule cleanupALL:
    input:
        raw = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/gene_counts_ALL.txt'

    output:
        cleaned = '/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_ALL.txt'

    shell:
        r'''
        cut -f 1,7-55 {input.raw} > /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/cut_ALL.txt

        tail -n +2 /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/cut_ALL.txt > /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/tail_ALL.txt

        awk -F'\t' 'BEGIN {{OFS=FS}} {{sub("/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/", "", $1); for (i=2; i<=NF; i++) sub("/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/", "", $i)}} NF>0' /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/tail_ALL.txt > /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_ALL.txt
        rm /home/workspace/jogrady/heQTL/work/RNA_seq/quantification/cut_ALL.txt
        '''

