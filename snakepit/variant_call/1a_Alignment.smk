import pandas as pd
import os
import subprocess

#configfile: 'config/config.yaml'
IDS_0 = [
    'SRR12609394', 'SRR12609405', 'SRR12609414', 'SRR12609425', 'SRR12609434',
    'SRR12609449', 'SRR12609475', 'SRR12609489', 'SRR12609495', 'SRR12609506',
    'SRR12609516', 'SRR12609526', 'SRR12609535', 'SRR12609546', 'SRR12609574',
    'SRR12609595', 'SRR12609608', 'SRR12609626', 'SRR12609641', 'SRR12609651',
    'SRR12609676', 'SRR12609687', 'SRR12609696', 'SRR12609770', 'SRR12609780',
    'SRR12609800', 'SRR12609808', 'SRR12609817', 'SRR12609855', 'SRR12609875',
    'SRR12609885', 'SRR12609897', 'SRR12609910', 'SRR12609922', 'SRR12609934',
    'SRR12609974', 'SRR12609982', 'SRR12609994', 'SRR12610003', 'SRR12610014',
    'SRR12610032', 'SRR12610053', 'SRR12610063', 'SRR12610071', 'SRR12610092',
    'SRR12610115', 'SRR12610123', 'SRR12610148']
IDS_1 = [
    'SRR12609772', 'SRR12609781', 'SRR12609801', 'SRR12609809', 'SRR12609818',
    'SRR12609395', 'SRR12609406', 'SRR12609415', 'SRR12609426', 'SRR12609435',
    'SRR12609856', 'SRR12609450', 'SRR12609866', 'SRR12609476', 'SRR12609876',
    'SRR12609886', 'SRR12609898', 'SRR12609911', 'SRR12609923', 'SRR12609935',
    'SRR12609496', 'SRR12609507', 'SRR12609517', 'SRR12609527', 'SRR12609536',
    'SRR12609547', 'SRR12609575', 'SRR12609975', 'SRR12609983', 'SRR12609995',
    'SRR12610004', 'SRR12610015', 'SRR12610054', 'SRR12609596', 'SRR12609609',
    'SRR12609627', 'SRR12609642', 'SRR12609677', 'SRR12610033', 'SRR12610064',
    'SRR12610072', 'SRR12609652', 'SRR12610093', 'SRR12610116', 'SRR12610124',
    'SRR12610149', 'SRR12609688', 'SRR12609697']
IDS_2 = [
    'SRR12609774', 'SRR12609783', 'SRR12609803', 'SRR12609812', 'SRR12609820',
    'SRR12609397', 'SRR12609409', 'SRR12609418', 'SRR12609429', 'SRR12609438',
    'SRR12609859', 'SRR12609453', 'SRR12609869', 'SRR12609480', 'SRR12609880',
    'SRR12609890', 'SRR12609901', 'SRR12609914', 'SRR12609926', 'SRR12609938',
    'SRR12609499', 'SRR12609509', 'SRR12609520', 'SRR12609530', 'SRR12609539',
    'SRR12609550', 'SRR12609578', 'SRR12609978', 'SRR12609986', 'SRR12609998',
    'SRR12610007', 'SRR12610018', 'SRR12610057', 'SRR12609599', 'SRR12609612',
    'SRR12609629', 'SRR12609645', 'SRR12609680', 'SRR12610036', 'SRR12610066',
    'SRR12610075', 'SRR12609655', 'SRR12610096', 'SRR12610119', 'SRR12610127',
    'SRR12610152', 'SRR12609691', 'SRR12609700']
IDS_3 = [
    'SRR12609775', 'SRR12609784', 'SRR12609804', 'SRR12609813', 'SRR12609821',
    'SRR12609399', 'SRR12609411', 'SRR12609420', 'SRR12609431', 'SRR12609439',
    'SRR12609861', 'SRR12609454', 'SRR12609871', 'SRR12609482', 'SRR12609882',
    'SRR12609892', 'SRR12609903', 'SRR12609916', 'SRR12609928', 'SRR12609940',
    'SRR12609501', 'SRR12609511', 'SRR12609522', 'SRR12609531', 'SRR12609541',
    'SRR12609551', 'SRR12609580', 'SRR12609980', 'SRR12609988', 'SRR12610000',
    'SRR12610009', 'SRR12610020', 'SRR12609587', 'SRR12609601', 'SRR12609614',
    'SRR12609634', 'SRR12609647', 'SRR12610058', 'SRR12610037', 'SRR12610068',
    'SRR12610077', 'SRR12609656', 'SRR12610098', 'SRR12610121', 'SRR12610129',
    'SRR12609682', 'SRR12609693', 'SRR12609702']
IDS_4 = [
    'SRR12609776', 'SRR12609787', 'SRR12609805', 'SRR12609814', 'SRR12609822',
    'SRR12609401', 'SRR12609412', 'SRR12609422', 'SRR12609433', 'SRR12609440',
    'SRR12609863', 'SRR12609457', 'SRR12609872', 'SRR12609484', 'SRR12609883',
    'SRR12609894', 'SRR12609905', 'SRR12609918', 'SRR12609930', 'SRR12609942',
    'SRR12609503', 'SRR12609513', 'SRR12609523', 'SRR12609533', 'SRR12609543',
    'SRR12609552', 'SRR12609581', 'SRR12609981', 'SRR12609990', 'SRR12610001',
    'SRR12610011', 'SRR12610021', 'SRR12609589', 'SRR12609603', 'SRR12609616',
    'SRR12609636', 'SRR12609649', 'SRR12610060', 'SRR12610038', 'SRR12610069',
    'SRR12610079', 'SRR12609658', 'SRR12610101', 'SRR12610122', 'SRR12610131',
    'SRR12609684', 'SRR12609694', 'SRR12609704']

rule all:
    input:
        'work/T0/fastqc/multiqc/multiqc_report.html',
        'work/T1/fastqc/multiqc/multiqc_report.html',
        'work/T2/fastqc/multiqc/multiqc_report.html',
        'work/T3/fastqc/multiqc/multiqc_report.html',
        'work/T4/fastqc/multiqc/multiqc_report.html',
        expand('work/T0/trimmed/{individual}_{N}_trimmed.fastq.gz', individual=config["fastqfiles_0"], N = (1,2)),
        expand('work/T1/trimmed/{individual}_{N}_trimmed.fastq.gz', individual=config["fastqfiles_1"], N = (1,2)),
        expand('work/T2/trimmed/{individual}_{N}_trimmed.fastq.gz', individual=config["fastqfiles_2"], N = (1,2)),
        expand('work/T3/trimmed/{individual}_{N}_trimmed.fastq.gz', individual=config["fastqfiles_3"], N = (1,2)),
        expand('work/T4/trimmed/{individual}_{N}_trimmed.fastq.gz', individual=config["fastqfiles_4"], N = (1,2)),
        'work/T0/fastqcTrim/multiqc/multiqc_report.html',
        'work/T1/fastqcTrim/multiqc/multiqc_report.html',
        'work/T2/fastqcTrim/multiqc/multiqc_report.html',
        'work/T3/fastqcTrim/multiqc/multiqc_report.html',
        'work/T4/fastqcTrim/multiqc/multiqc_report.html',
        expand('work/T0/align/{individual}.bam.bai', individual=config["fastqfiles_0"]),
        expand('work/T1/align/{individual}.bam.bai', individual=config["fastqfiles_1"]),
        expand('work/T2/align/{individual}.bam.bai', individual=config["fastqfiles_2"]),
        expand('work/T3/align/{individual}.bam.bai', individual=config["fastqfiles_3"]),
        expand('work/T4/align/{individual}.bam.bai', individual=config["fastqfiles_4"])

rule fastqc_0:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_0"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('work/T0/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/T0/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/T0/fastqc/'


rule fastqc_1:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_1"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('work/T1/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/T1/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/T1/fastqc/'


rule fastqc_2:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_2"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('work/T2/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/T2/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/T2/fastqc/'


rule fastqc_3:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_3"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('work/T3/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/T3/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/T3/fastqc/'

rule fastqc_4:
    input:
        reads = lambda wildcards: expand(f'{config["fastqfiles_4"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('work/T4/fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('work/T4/fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/T4/fastqc/'


rule multiqc_0:
    input:
        expand(['work/T0/fastqc/{srr}_1_fastqc.zip', 'work/T0/fastqc/{srr}_2_fastqc.zip'], srr = IDS_0)
    output:
        report='work/T0/fastqc/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T0/fastqc/multiqc/
        """

rule multiqc_1:
    input:
        expand(['work/T1/fastqc/{srr}_1_fastqc.zip', 'work/T1/fastqc/{srr}_2_fastqc.zip'], srr = IDS_1)
    output:
        report='work/T1/fastqc/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T1/fastqc/multiqc/
        """
rule multiqc_2:
    input:
        expand(['work/T2/fastqc/{srr}_1_fastqc.zip', 'work/T2/fastqc/{srr}_2_fastqc.zip'], srr = IDS_2)
    output:
        report='work/T2/fastqc/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T2/fastqc/multiqc/
        """
rule multiqc_3:
    input:
        expand(['work/T3/fastqc/{srr}_1_fastqc.zip', 'work/T3/fastqc/{srr}_2_fastqc.zip'], srr = IDS_3)
    output:
        report='work/T3/fastqc/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T3/fastqc/multiqc/
        """
rule multiqc_4:
    input:
        expand(['work/T4/fastqc/{srr}_1_fastqc.zip', 'work/T4/fastqc/{srr}_2_fastqc.zip'], srr = IDS_4)
    output:
        report='work/T4/fastqc/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T4/fastqc/multiqc/
        """



rule trimming_0:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_0"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/T0/trimmed/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/T0/trimmed/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'


rule trimming_1:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_1"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/T1/trimmed/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/T1/trimmed/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'

rule trimming_2:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_2"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/T2/trimmed/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/T2/trimmed/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'


rule trimming_3:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_3"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/T3/trimmed/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/T3/trimmed/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'

rule trimming_4:
    input:
        untrimmed_reads = lambda wildcards: expand(f'{config["fastqfiles_4"][wildcards.individual]}_{{N}}.fastq.gz', N = (1,2)),
        adapters = config["adapter"]["illumina"]
    output:
        trimmed_reads=expand('work/T4/trimmed/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/T4/trimmed/{{individual}}_{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 MINLEN:36 TRAILING:30'



rule fastqcTrimmed_0:
    input:
        Tread = lambda wildcards: expand(f'work/T0/trimmed/{config["IDS_0"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/T0/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/T0/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/T0/fastqcTrim/'
        
rule fastqcTrimmed_1:
    input:
        Tread = lambda wildcards: expand(f'work/T1/trimmed/{config["IDS_1"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/T1/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/T1/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/T1/fastqcTrim/'

rule fastqcTrimmed_2:
    input:
        Tread = lambda wildcards: expand(f'work/T2/trimmed/{config["IDS_2"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/T2/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/T2/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/T2/fastqcTrim/'
        


rule fastqcTrimmed_3:
    input:
        Tread = lambda wildcards: expand(f'work/T3/trimmed/{config["IDS_3"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/T3/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/T3/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/T3/fastqcTrim/'
        
        
        
rule fastqcTrimmed_4:
    input:
        Tread = lambda wildcards: expand(f'work/T4/trimmed/{config["IDS_4"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2))
    output:
        reads_trimmed = expand('work/T4/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.zip', N=(1,2)),
        html_trimmed = expand('work/T4/fastqcTrim/{{srr}}_{N}_trimmed_fastqc.html', N=(1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.Tread[0]} {input.Tread[1]} -t {threads} -o work/T4/fastqcTrim/'

rule multiqcTrimmed_0:
    input:
        expand(['work/T0/fastqcTrim/{srr}_1_trimmed_fastqc.zip', 'work/T0/fastqcTrim/{srr}_2_trimmed_fastqc.zip'], srr = IDS_0)
    output:
        report='work/T0/fastqcTrim/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T0/fastqcTrim/multiqc/
        """

rule multiqcTrimmed_1:
    input:
        expand(['work/T1/fastqcTrim/{srr}_1_trimmed_fastqc.zip', 'work/T1/fastqcTrim/{srr}_2_trimmed_fastqc.zip'], srr = IDS_1)
    output:
        report='work/T1/fastqcTrim/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T1/fastqcTrim/multiqc/
        """
rule multiqcTrimmed_2:
    input:
        expand(['work/T2/fastqcTrim/{srr}_1_trimmed_fastqc.zip', 'work/T2/fastqcTrim/{srr}_2_trimmed_fastqc.zip'], srr = IDS_2)
    output:
        report='work/T2/fastqcTrim/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T2/fastqcTrim/multiqc/
        """
rule multiqcTrimmed_3:
    input:
        expand(['work/T3/fastqcTrim/{srr}_1_trimmed_fastqc.zip', 'work/T3/fastqcTrim/{srr}_2_trimmed_fastqc.zip'], srr = IDS_3)
    output:
        report='work/T3/fastqcTrim/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T3/fastqcTrim/multiqc/
        """
rule multiqcTrimmed_4:
    input:
        expand(['work/T4/fastqcTrim/{srr}_1_trimmed_fastqc.zip', 'work/T4/fastqcTrim/{srr}_2_trimmed_fastqc.zip'], srr = IDS_4)
    output:
        report='work/T4/fastqcTrim/multiqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/T4/fastqcTrim/multiqc/
        """


rule GenomeBuild:
    input:
        fasta = config["starinputs"]["reference"],
        annotation = config["starinputs"]["annot"]
    output:
        directory('work/STAR_DIR')

    conda:
        'envs/environment.yaml'
    threads: 20
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.annotation} '
        '--sjdbOverhang 99'


rule Alignment_0:
    input:
        genome = 'work/STAR_DIR',
        reads = lambda wildcards: expand(f'work/T0/trimmed/{config["IDS_0"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_0"][wildcards.srr]}'
    output:
        aligned = 'work/T0/align/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/T0/align/{srr}_Log.final.out',
        interlog = 'work/T0/align/{srr}_Log.progress.out',
        initiallog = 'work/T0/align/{srr}_Log.out'
    conda:
        'envs/environment.yaml'

    threads: 20

    shell:
        'STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/T0/align/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 10000000000'


rule Alignment_1:
    input:
        genome = 'work/STAR_DIR',
        reads = lambda wildcards: expand(f'work/T1/trimmed/{config["IDS_1"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_1"][wildcards.srr]}'
    output:
        aligned = 'work/T1/align/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/T1/align/{srr}_Log.final.out',
        interlog = 'work/T1/align/{srr}_Log.progress.out',
        initiallog = 'work/T1/align/{srr}_Log.out'
    conda:
        'envs/environment.yaml'

    threads: 20


    shell:
        'STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/T1/align/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 50000000000'
        
        
rule Alignment_2:
    input:
        genome = 'work/STAR_DIR',
        reads = lambda wildcards: expand(f'work/T2/trimmed/{config["IDS_2"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_2"][wildcards.srr]}'
    output:
        aligned = 'work/T2/align/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/T2/align/{srr}_Log.final.out',
        interlog = 'work/T2/align/{srr}_Log.progress.out',
        initiallog = 'work/T2/align/{srr}_Log.out'

    conda:
        'envs/environment.yaml'

    threads: 20


    shell:
        'STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/T2/align/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 50000000000'
  
  
  
rule Alignment_3:
    input:
        genome = 'work/STAR_DIR',
        reads = lambda wildcards: expand(f'work/T3/trimmed/{config["IDS_3"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_3"][wildcards.srr]}'
    output:
        aligned = 'work/T3/align/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/T3/align/{srr}_Log.final.out',
        interlog = 'work/T3/align/{srr}_Log.progress.out',
        initiallog = 'work/T3/align/{srr}_Log.out'
    conda:
        'envs/environment.yaml'

    threads: 20


    shell:
        'STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/T3/align/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 50000000000'
        

rule Alignment_4:
    input:
        genome = 'work/STAR_DIR',
        reads = lambda wildcards: expand(f'work/T4/trimmed/{config["IDS_4"][wildcards.srr]}_{{N}}_trimmed.fastq.gz', N=(1,2)),

    params:
        prefix = lambda wildcards: f'{config["IDS_4"][wildcards.srr]}'
    output:
        aligned = 'work/T4/align/{srr}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/T4/align/{srr}_Log.final.out',
        interlog = 'work/T4/align/{srr}_Log.progress.out',
        initiallog = 'work/T4/align/{srr}_Log.out'
    conda:
        'envs/environment.yaml'

    threads: 20


    shell:
        'STAR --genomeLoad NoSharedMemory '
        '--genomeDir {input.genome} --runThreadN {threads} '
        '--readFilesIn {input.reads[0]} {input.reads[1]} '
        '--readFilesCommand gunzip -c '
        '--outFileNamePrefix work/T4/align/{params.prefix}_ '
        '--outSAMtype BAM SortedByCoordinate '
        '--limitBAMsortRAM 50000000000'


rule renamebam_0:
    input:
        bam = expand('work/T0/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_0)
    output:
        simple_index = 'work/T0/align/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
        'samtools index {input.bam}'

rule renamebam_1:
    input:
        bam = expand('work/T1/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_1)
    output:
        simple_index = 'work/T1/align/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
        'samtools index {input.bam}'
rule renamebam_2:
    input:
        bam = expand('work/T2/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_2)
    output:
        simple_index = 'work/T2/align/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
         'samtools index {input.bam}'
rule renamebam_3:
    input:
        bam = expand('work/T3/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_3)
    output:
        simple_index = 'work/T3/align/{srr}_Aligned.sortedByCoord.out.bam.bai'
    threads: 4

    shell:
         'samtools index {input.bam}'
rule renamebam_4:
    input:
        bam = expand('work/T4/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_4)
    output:
        simple_index = 'work/T4/align/{srr}_Aligned.sortedByCoord.out.bam.bai',
    threads: 4

    shell:
         'samtools index {input.bam}'


rule clean_0:
    input:
        bam = expand('work/T0/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_0),
        bai = expand('work/T0/align/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = IDS_0)
    output:
        simple_bam = 'work/T0/align/{srr}.bam',
        simple_bai = 'work/T0/align/{srr}.bam.bai'
    shell:
         'rename "_Aligned.sortedByCoord.out" "" {input.bam} && '
         'rename "_Aligned.sortedByCoord.out" "" {input.bai}'

rule clean_1:
    input:
        bam = expand('work/T1/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_1),
        bai = expand('work/T1/align/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = IDS_1)
    output:
        simple_bam = 'work/T1/align/{srr}.bam',
        simple_bai = 'work/T1/align/{srr}.bam.bai'
    shell:
         'rename "_Aligned.sortedByCoord.out" "" {input.bam} && '
         'rename "_Aligned.sortedByCoord.out" "" {input.bai}'

rule clean_2:
    input:
        bam = expand('work/T2/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_2),
        bai = expand('work/T2/align/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = IDS_2)
    output:
        simple_bam = 'work/T2/align/{srr}.bam',
        simple_bai = 'work/T2/align/{srr}.bam.bai'
    shell:
         'rename "_Aligned.sortedByCoord.out" "" {input.bam} && '
         'rename "_Aligned.sortedByCoord.out" "" {input.bai}'

rule clean_3:
    input:
        bam = expand('work/T3/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_3),
        bai = expand('work/T3/align/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = IDS_3)
    output:
        simple_bam = 'work/T3/align/{srr}.bam',
        simple_bai = 'work/T3/align/{srr}.bam.bai'
    shell:
         'rename "_Aligned.sortedByCoord.out" "" {input.bam} && '
         'rename "_Aligned.sortedByCoord.out" "" {input.bai}'

rule clean_4:
    input:
        bam = expand('work/T4/align/{{srr}}_Aligned.sortedByCoord.out.bam', srr = IDS_4),
        bai = expand('work/T4/align/{{srr}}_Aligned.sortedByCoord.out.bam.bai', srr = IDS_4)
    output:
        simple_bam = 'work/T4/align/{srr}.bam',
        simple_bai = 'work/T4/align/{srr}.bam.bai'
    shell:
         'rename "_Aligned.sortedByCoord.out" "" {input.bam} && '
         'rename "_Aligned.sortedByCoord.out" "" {input.bai}'

