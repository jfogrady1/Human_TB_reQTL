#rule all:
    #input:
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage.mosdepth.region.dist.txt', id=config["merged_patient_IDS"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/coding_regions.html'),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage_WG.mosdepth.global.dist.txt', id=config["merged_patient_IDS"]),
        #expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/WG_regions.html')


rule CDS:
    input:
        annotation = '/home/workspace/jogrady/heQTL/data/GIAB/gencode.v43.basic.annotation.gff3.gz'

    output:
        cds = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/CDS.bed'

    shell:
        '''
        gzip -dc {input.annotation} | awk -v OFS='\t' '$3 == "CDS"  && $4 < $5 {{ print $1, $4, $5, "CDS" }}' | awk '!dup[$0]++' > {output.cds}
        '''

rule mosdepth_coverage:
    input:
        bam = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}.bam',
        cds = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/CDS.bed'
    output:
        coverage = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage.mosdepth.region.dist.txt'
    threads: 8
    shell:
        '''
        /home/workspace/jogrady/heQTL/software/mosdepth-0.3.9 --by {input.cds} --threads {threads} /home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{wildcards.id}_coverage {input.bam}
        '''

rule plotting:
    input:
        cov = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage.mosdepth.region.dist.txt', id=config["merged_patient_IDS"]),
        script = '/home/workspace/jogrady/heQTL/scripts/mosdepth_plot.py'
    output:
        html = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/coding_regions.html'

    shell:
        '''
        python {input.script} -o {output.html} {input.cov}
        '''     

rule mosdepth_coverage_WG:
    input:
        bam = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}.bam',
    output:
        coverage = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage_WG.mosdepth.global.dist.txt'
    threads: 10
    shell:
        '''
        /home/workspace/jogrady/heQTL/software/mosdepth-0.3.9 --threads {threads} /home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{wildcards.id}_coverage_WG {input.bam}
        '''


rule plotting_WG:
    input:
        cov = expand('/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/{id}_coverage_WG.mosdepth.global.dist.txt', id=config["merged_patient_IDS"]),
        script = '/home/workspace/jogrady/heQTL/scripts/mosdepth_plot.py'
    output:
        html = '/home/workspace/jogrady/heQTL/work/RNA_seq/align/mergedfinal/WG_regions.html'

    shell:
        '''
        python {input.script} -o {output.html} {input.cov}
        '''

