import shutil
import glob

P_3IDS = ['SRR12609780','SRR12609784',	'SRR12609787']

# No timepoint 0 for sample P_157
P_157_IDs = ['SRR12609496',	'SRR12609499',	'SRR12609501',	'SRR12609503']

cohort_samples0 = config["SRR_merge_0"] 
cohort_samples1 = config["SRR_merge_1"] 
cohort_samples2 = config["SRR_merge_2"] 
cohort_samples3 = config["SRR_merge_3"]
cohort_samples4 = config["SRR_merge_4"] 
cohort_ids = config["patient_IDS"]
rule all:
    input:
       ['/cluster/work/pausch/jogrady/heQTL/work/mergedbam/{}_{}_{}_{}_{}/{}.bam'.format(u,v,w,x,y,z) for u,v,w,x,y,z in zip(cohort_samples0, cohort_samples1, cohort_samples2, cohort_samples3, cohort_samples4, cohort_ids)],
       '/cluster/work/pausch/jogrady/heQTL/work/mergedbam/SRR12609780_SRR12609783_SRR12609784_SRR12609787/P_3.bam',
       '/cluster/work/pausch/jogrady/heQTL/work/mergedbam/SRR12609496_SRR12609499_SRR12609501_SRR12609503/P_157.bam'
rule mergebam:
    input:
        T0 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T0/align/{config["SRR_merge_0"][wildcards.srr0]}.bam'),
        T1 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T1/align/{config["SRR_merge_1"][wildcards.srr1]}.bam'),
        T2 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T2/align/{config["SRR_merge_2"][wildcards.srr2]}.bam'), 
        T3 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T3/align/{config["SRR_merge_3"][wildcards.srr3]}.bam'), 
        T4 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T4/align/{config["SRR_merge_4"][wildcards.srr4]}.bam'),
    output:
        merged_temp = '/cluster/work/pausch/jogrady/heQTL/work/mergedbam/{srr0}_{srr1}_{srr2}_{srr3}_{srr4}/{id}.bam'
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = "8:00",
        disk_scratch = 30
    shell:
      '''
         samtools merge --threads {threads} -o {output.merged_temp} {input.T0} {input.T1} {input.T2} {input.T3} {input.T4}

      '''
rule merge3:
    input:
        T0 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T0/align/SRR12609780.bam'),
        T2 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T2/align/SRR12609783.bam'),
        T3 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T3/align/SRR12609784.bam'),
        T4 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T4/align/SRR12609787.bam')

    output:
        merged_temp = '/cluster/work/pausch/jogrady/heQTL/work/mergedbam/SRR12609780_SRR12609783_SRR12609784_SRR12609787/P_3.bam'
    shell:
        '''
         samtools merge --threads {threads} -o {output.merged_temp} {input.T0} {input.T2} {input.T3} {input.T4}

        '''

rule merge157:
    input:
        T1 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T1/align/SRR12609496.bam'),
        T2 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T2/align/SRR12609499.bam'),
        T3 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T3/align/SRR12609501.bam'),
        T4 = lambda wildcards: expand(f'/cluster/work/pausch/jogrady/heQTL/work/T4/align/SRR12609503.bam')

    output:
        merged_temp = '/cluster/work/pausch/jogrady/heQTL/work/mergedbam/SRR12609496_SRR12609499_SRR12609501_SRR12609503/P_157.bam'
    shell:
        '''
         samtools merge --threads {threads} -o {output.merged_temp} {input.T1} {input.T2} {input.T3} {input.T4}
        '''
