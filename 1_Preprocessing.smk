import pandas as pd
import os
import subprocess

configfile: "config/Alignment.yaml"
# File to preprocess RNA-seq data and genotype data
#timepoint = [0,1,2,3,4]
reads = ["1","2"]
autosomes = [str(i) for i in range(1,22)] # bovine autosomes

# Get sample names into a list
#T0
samples_information = pd.read_csv("./data/T0/accession_file.txt", sep='\t', index_col=False, header = None)
sample_names_0 = list(samples_information[0])

#T1
samples_information = pd.read_csv("./data/T1/accession_file.txt", sep='\t', index_col=False, header = None)
sample_names_1 = list(samples_information[0])

#T2
samples_information = pd.read_csv("./data/T2/accession_file.txt", sep='\t', index_col=False, header = None)
sample_names_2 = list(samples_information[0])

#T3
samples_information = pd.read_csv("./data/T3/accession_file.txt", sep='\t', index_col=False, header = None)
sample_names_3 = list(samples_information[0])

#T4
samples_information = pd.read_csv("./data/T4/accession_file.txt", sep='\t', index_col=False, header = None)
sample_names_4 = list(samples_information[0])


rule all:
    input:
        expand('data/T4/{SRR_4}_{n}.fastq.gz', n = reads, SRR_4 = config["IDS_4"]),
        expand('work/RNA_seq/QC/fastqc/{runs}_{n}_fastqc.zip', runs = config["ALL_samples"], n = reads),
        'work/RNA_seq/trimmed/fastqc/multiqc_report.html',
        expand('work/RNA_seq/trimmed/T0/{individual}_{N}_unpaired.fastq.gz', N = (1,2), individual = config["IDS_0"]),
        expand('work/RNA_seq/trimmed/T1/{individual}_{N}_unpaired.fastq.gz', N = (1,2), individual = config["IDS_1"]),
        expand('work/RNA_seq/trimmed/T2/{individual}_{N}_unpaired.fastq.gz', N = (1,2), individual = config["IDS_2"]),
        expand('work/RNA_seq/trimmed/T3/{individual}_{N}_unpaired.fastq.gz', N = (1,2), individual = config["IDS_3"]),
        expand('work/RNA_seq/trimmed/T4/{individual}_{N}_unpaired.fastq.gz', N = (1,2), individual = config["IDS_4"]),
        directory('work/RNA_seq/STAR_DIR'),
        expand('work/RNA_seq/align/T0/{srr}_Log.final.out', srr = config["IDS_0"]),
        expand('work/RNA_seq/align/T1/{srr}_Log.final.out', srr = config["IDS_1"]),
        expand('work/RNA_seq/align/T2/{srr}_Log.final.out', srr = config["IDS_2"]),
        expand('work/RNA_seq/align/T3/{srr}_Log.final.out', srr = config["IDS_3"]),
        expand('work/RNA_seq/align/T4/{srr}_Log.final.out', srr = config["IDS_4"]),
        expand('work/RNA_seq/align/T0/{srr}.bam', srr = config["IDS_0"]),
        expand('work/RNA_seq/align/T1/{srr}.bam', srr = config["IDS_1"]),
        expand('work/RNA_seq/align/T2/{srr}.bam', srr = config["IDS_2"]),
        expand('work/RNA_seq/align/T3/{srr}.bam', srr = config["IDS_3"]),
        expand('work/RNA_seq/align/T4/{srr}.bam', srr = config["IDS_4"])
        
              



rule download:
    output:
        srr_0 = expand("data/T0/{SRR_0}_{n}.fastq.gz", SRR_0 = config["IDS_0"], n = reads),
        srr_1 = expand("data/T1/{SRR_1}_{n}.fastq.gz", SRR_1 = config["IDS_1"], n = reads),
        srr_2 = expand("data/T2/{SRR_2}_{n}.fastq.gz", SRR_2 = config["IDS_2"], n = reads),
        srr_3 = expand("data/T3/{SRR_3}_{n}.fastq.gz", SRR_3 = config["IDS_3"], n = reads),
        srr_4 = expand("data/T4/{SRR_4}_{n}.fastq.gz", SRR_4 = config["IDS_4"], n = reads)

    conda:
        "fastq-dl"
    
    params:
        samples_0 = sample_names_0,
        samples_1 = sample_names_1,
        samples_2 = sample_names_2,
        samples_3 = sample_names_3,
        samples_4 = sample_names_4

    shell:
        '''
        for i in {params.samples_0}
        do
        fastq-dl --cpus 10 -a $i --outdir data/T0/; 
        done
        


        for i in {params.samples_1}
        do
        fastq-dl --cpus 10 -a $i --outdir data/T1/; 
        done
        


        for i in {params.samples_2}
        do
        fastq-dl --cpus 10 -a $i --outdir data/T2/; 
        done
        


        for i in {params.samples_3}
        do
        fastq-dl --cpus 10 -a $i --outdir data/T3/; 
        done


        for i in {params.samples_4}
        do
        fastq-dl --cpus 10 -a $i --outdir data/T4/; 
        done
        '''

rule fastqc0:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_0"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        zip_file = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads: 20
    shell:
        '''
        fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/QC/fastqc/
        '''

rule fastqc1:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_1"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        zip_file = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads: 20
    shell:
        '''
        fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/QC/fastqc/
        '''


rule fastqc2:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_2"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        zip_file = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads: 20
    shell:
        '''
        fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/QC/fastqc/
        '''


rule fastqc3:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_3"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        zip_file = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads: 20
    shell:
        '''
        fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/QC/fastqc/
        '''

rule fastqc4:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_4"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        zip_file = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/QC/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads: 20
    shell:
        '''
        fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/QC/fastqc/
        '''

rule multiqc:
    input:
        reads = expand('work/RNA_seq/QC/fastqc/{individual}_{N}_fastqc.zip', individual = config['ALL_samples'], N = (1,2))
    output:
        report='work/RNA_seq/QC/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/QC/fastqc/
        """

rule trimming_0:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_0"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/RNA_seq/trimmed/T0/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/trimmed/T0/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    singularity:
        "docker://staphb/trimmomatic:latest"
    resources:
        mem_mb = 4000
    shell:
        '/Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'

rule trimming_1:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_1"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/RNA_seq/trimmed/T1/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/trimmed/T1/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    singularity:
        "docker://staphb/trimmomatic:latest"
    resources:
        mem_mb = 4000
    shell:
        '/Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'




rule trimming_2:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_2"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/RNA_seq/trimmed/T2/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/trimmed/T2/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    singularity:
        "docker://staphb/trimmomatic:latest"
    resources:
        mem_mb = 4000
    shell:
        '/Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'


rule trimming_3:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_3"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/RNA_seq/trimmed/T3/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/trimmed/T3/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    singularity:
        "docker://staphb/trimmomatic:latest"
    resources:
        mem_mb = 4000
    shell:
        '/Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'

rule trimming_4:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_4"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/RNA_seq/trimmed/T4/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/trimmed/T4/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    singularity:
        "docker://staphb/trimmomatic:latest"
    resources:
        mem_mb = 4000
    shell:
        '/Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'


rule fastqcTrimmed_0:
    input:
        Tread = lambda wildcards: expand(f'work/RNA_seq/trimmed/T0/{config["IDS_0"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/RNA_seq/trimmed/fastqc/'

rule fastqcTrimmed_1:
    input:
        Tread = lambda wildcards: expand(f'work/RNA_seq/trimmed/T1/{config["IDS_1"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/RNA_seq/trimmed/fastqc/'

rule fastqcTrimmed_2:
    input:
        Tread = lambda wildcards: expand(f'work/RNA_seq/trimmed/T2/{config["IDS_2"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/RNA_seq/trimmed/fastqc/'

rule fastqcTrimmed_3:
    input:
        Tread = lambda wildcards: expand(f'work/RNA_seq/trimmed/T3/{config["IDS_3"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/RNA_seq/trimmed/fastqc/'

rule fastqcTrimmed_4:
    input:
        Tread = lambda wildcards: expand(f'work/RNA_seq/trimmed/T4/{config["IDS_4"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/RNA_seq/trimmed/fastqc/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/RNA_seq/trimmed/fastqc/'


rule multiqcTrimmed:
    input:
        expand(['work/RNA_seq/trimmed/fastqc/{srr}_1_trimmed_fastqc.zip', 'work/RNA_seq/trimmed/fastqc/{srr}_2_trimmed_fastqc.zip'], srr = config["ALL_samples"])
    output:
        report='work/RNA_seq/trimmed/fastqc/multiqc_report.html'
        
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/trimmed/fastqc/
        """

rule GenomeBuild:
    input:
        fasta = config["starinputs"]["reference"],
        annotation = config["starinputs"]["annot"]
    output:
        star_dir = directory('work/RNA_seq/STAR_DIR')
    threads: 20
    shell:
        'mkdir {output} && '
        'software/STAR-2.7.10b/source/STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.annotation} '
        '--sjdbOverhang 99'


rule Alignment_0:
    input:
        genome = rules.GenomeBuild.output.star_dir,
        reads = lambda wildcards: expand(f'work/RNA_seq/trimmed/T0/{config["IDS_0"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_0"][wildcards.srr]}'
    output:
        aligned = 'work/RNA_seq/align/T0/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/align/T0/{srr}_Log.final.out',
        interlog = 'work/RNA_seq/align/T0/{srr}_Log.progress.out',
        initiallog = 'work/RNA_seq/align/T0/{srr}_Log.out'
    threads: 20

    shell:
        'software/STAR-2.7.10b/source/STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/RNA_seq/align/T0/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 10000000000'


rule Alignment_1:
    input:
        genome = rules.GenomeBuild.output.star_dir,
        reads = lambda wildcards: expand(f'work/RNA_seq/trimmed/T1/{config["IDS_1"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_1"][wildcards.srr]}'
    output:
        aligned = 'work/RNA_seq/align/T1/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/align/T1/{srr}_Log.final.out',
        interlog = 'work/RNA_seq/align/T1/{srr}_Log.progress.out',
        initiallog = 'work/RNA_seq/align/T1/{srr}_Log.out'
    threads: 20

    shell:
        'software/STAR-2.7.10b/source/STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/RNA_seq/align/T1/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 10000000000'


rule Alignment_2:
    input:
        genome = rules.GenomeBuild.output.star_dir,
        reads = lambda wildcards: expand(f'work/RNA_seq/trimmed/T2/{config["IDS_2"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_2"][wildcards.srr]}'
    output:
        aligned = 'work/RNA_seq/align/T2/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/align/T2/{srr}_Log.final.out',
        interlog = 'work/RNA_seq/align/T2/{srr}_Log.progress.out',
        initiallog = 'work/RNA_seq/align/T2/{srr}_Log.out'
    threads: 20

    shell:
        'software/STAR-2.7.10b/source/STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/RNA_seq/align/T2/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 10000000000'



rule Alignment_3:
    input:
        genome = rules.GenomeBuild.output.star_dir,
        reads = lambda wildcards: expand(f'work/RNA_seq/trimmed/T3/{config["IDS_3"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_3"][wildcards.srr]}'
    output:
        aligned = 'work/RNA_seq/align/T3/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/align/T3/{srr}_Log.final.out',
        interlog = 'work/RNA_seq/align/T3/{srr}_Log.progress.out',
        initiallog = 'work/RNA_seq/align/T3/{srr}_Log.out'
    threads: 20

    shell:
        'software/STAR-2.7.10b/source/STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/RNA_seq/align/T3/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 10000000000'




rule Alignment_4:
    input:
        genome = rules.GenomeBuild.output.star_dir,
        reads = lambda wildcards: expand(f'work/RNA_seq/trimmed/T4/{config["IDS_4"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_4"][wildcards.srr]}'
    output:
        aligned = 'work/RNA_seq/align/T4/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/align/T4/{srr}_Log.final.out',
        interlog = 'work/RNA_seq/align/T4/{srr}_Log.progress.out',
        initiallog = 'work/RNA_seq/align/T4/{srr}_Log.out'
    threads: 20

    shell:
        'software/STAR-2.7.10b/source/STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/RNA_seq/align/T4/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 10000000000'

rule index_bam_0:
    input:
        bam = expand('work/RNA_seq/align/T0/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_0"])
    output:
        simple_index = 'work/RNA_seq/align/T0/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
        'software/samtools-1.15.1/samtools index {input.bam}'

rule index_bam_1:
    input:
        bam = expand('work/RNA_seq/align/T1/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_1"])
    output:
        simple_index = 'work/RNA_seq/align/T1/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
        'software/samtools-1.15.1/samtools index {input.bam}'
rule index_bam_2:
    input:
        bam = expand('work/RNA_seq/align/T2/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_2"])
    output:
        simple_index = 'work/RNA_seq/align/T2/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
         'software/samtools-1.15.1/samtools index {input.bam}'
rule index_bam_3:
    input:
        bam = expand('work/RNA_seq/align/T3/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_3"])
    output:
        simple_index = 'work/RNA_seq/align/T3/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
         'software/samtools-1.15.1/samtools index {input.bam}'
rule index_bam_4:
    input:
        bam = expand('work/RNA_seq/align/T4/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_4"])
    output:
        simple_index = 'work/RNA_seq/align/T4/{srr}_Aligned.sortedByCoord.out.bam.bai',
    threads: 4

    shell:
         'software/samtools-1.15.1/samtools index {input.bam}'


rule rename_0:
    input:
        bam = expand('work/RNA_seq/align/T0/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_0"]),
        bai = expand('work/RNA_seq/align/T0/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = config["IDS_0"])
    output:
        simple_bam = 'work/RNA_seq/align/T0/{srr}.bam',
        simple_bai = 'work/RNA_seq/align/T0/{srr}.bam.bai'
    run:
        import os

        # Rename the BAM file
        os.rename(str(input.bam[0]), str(output.simple_bam))

        # Rename the BAI file
        os.rename(str(input.bai[0]), str(output.simple_bai))

rule rename_1:
    input:
        bam = expand('work/RNA_seq/align/T1/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_1"]),
        bai = expand('work/RNA_seq/align/T1/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = config["IDS_1"])
    output:
        simple_bam = 'work/RNA_seq/align/T1/{srr}.bam',
        simple_bai = 'work/RNA_seq/align/T1/{srr}.bam.bai'
    run:
        import os

        # Rename the BAM file
        os.rename(str(input.bam[0]), str(output.simple_bam))

        # Rename the BAI file
        os.rename(str(input.bai[0]), str(output.simple_bai))

rule rename_2:
    input:
        bam = expand('work/RNA_seq/align/T2/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_2"]),
        bai = expand('work/RNA_seq/align/T2/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = config["IDS_2"])
    output:
        simple_bam = 'work/RNA_seq/align/T2/{srr}.bam',
        simple_bai = 'work/RNA_seq/align/T2/{srr}.bam.bai'
    run:
        import os

        # Rename the BAM file
        os.rename(str(input.bam[0]), str(output.simple_bam))

        # Rename the BAI file
        os.rename(str(input.bai[0]), str(output.simple_bai))

rule rename_3:
    input:
        bam = expand('work/RNA_seq/align/T3/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_3"]),
        bai = expand('work/RNA_seq/align/T3/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = config["IDS_3"])
    output:
        simple_bam = 'work/RNA_seq/align/T3/{srr}.bam',
        simple_bai = 'work/RNA_seq/align/T3/{srr}.bam.bai'
    run:
        import os

        # Rename the BAM file
        os.rename(str(input.bam[0]), str(output.simple_bam))

        # Rename the BAI file
        os.rename(str(input.bai[0]), str(output.simple_bai))

rule rename_4:
    input:
        bam = expand('work/RNA_seq/align/T4/{{srr}}_Aligned.sortedByCoord.out.bam', srr = config["IDS_4"]),
        bai = expand('work/RNA_seq/align/T4/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = config["IDS_4"])
    output:
        simple_bam = 'work/RNA_seq/align/T4/{srr}.bam',
        simple_bai = 'work/RNA_seq/align/T4/{srr}.bam.bai'
    run:
        import os

        # Rename the BAM file
        os.rename(str(input.bam[0]), str(output.simple_bam))

        # Rename the BAI file
        os.rename(str(input.bai[0]), str(output.simple_bai))
