import shutil
import glob

configfile: "config/Alignment.yaml"

P_3IDS = ['SRR12609780','SRR12609784',	'SRR12609787']

# No timepoint 0 for sample P_157
P_157_IDs = ['SRR12609496',	'SRR12609499',	'SRR12609501',	'SRR12609503']

cohort_samples0 = config["SRR_merge_0"] 
cohort_samples1 = config["SRR_merge_1"] 
cohort_samples2 = config["SRR_merge_2"] 
cohort_samples3 = config["SRR_merge_3"]
cohort_samples4 = config["SRR_merge_4"] 
cohort_ids = config["patient_IDS"]
#rule all:
 #   input:
  #      'work/RNA_seq/align/mergedfinal/P_1.bam',
   #     'work/RNA_seq/align/mergedfinal/P_75.bam',
    #    'work/RNA_seq/align/mergedfinal/P_76.bam',
     #   'work/RNA_seq/align/mergedfinal/P_83.bam',
      #  'work/RNA_seq/align/mergedfinal/P_84.bam',
       # 'work/RNA_seq/align/mergedfinal/P_85.bam',
        #'work/RNA_seq/align/mergedfinal/P_92.bam',
        #'work/RNA_seq/align/mergedfinal/P_105.bam',
        #'work/RNA_seq/align/mergedfinal/P_99.bam',
        #'work/RNA_seq/align/mergedfinal/P_157.bam',
        #'work/RNA_seq/align/mergedfinal/P_175.bam',
        #'work/RNA_seq/align/mergedfinal/P_176.bam',
        #'work/RNA_seq/align/mergedfinal/P_201.bam',
        #'work/RNA_seq/align/mergedfinal/P_203.bam',
        #'work/RNA_seq/align/mergedfinal/P_221.bam',
        #'work/RNA_seq/align/mergedfinal/P_227.bam',
        #'work/RNA_seq/align/mergedfinal/P_268.bam',
        #'work/RNA_seq/align/mergedfinal/P_278.bam',
        #'work/RNA_seq/align/mergedfinal/P_294.bam',
        #'work/RNA_seq/align/mergedfinal/P_297.bam',
        #'work/RNA_seq/align/mergedfinal/P_367.bam',
        #'work/RNA_seq/align/mergedfinal/P_302.bam',
        #'work/RNA_seq/align/mergedfinal/P_432.bam',
        #'work/RNA_seq/align/mergedfinal/P_461.bam',
        #'work/RNA_seq/align/mergedfinal/P_1.bam',
        #'work/RNA_seq/align/mergedfinal/P_3.bam',
        #'work/RNA_seq/align/mergedfinal/P_13.bam',
        #'work/RNA_seq/align/mergedfinal/P_29.bam',
        #'work/RNA_seq/align/mergedfinal/P_36.bam',
        #'work/RNA_seq/align/mergedfinal/P_87.bam',
        #'work/RNA_seq/align/mergedfinal/P_117.bam',
        #'work/RNA_seq/align/mergedfinal/P_124.bam',
        #'work/RNA_seq/align/mergedfinal/P_127.bam',
        #'work/RNA_seq/align/mergedfinal/P_136.bam',
        #'work/RNA_seq/align/mergedfinal/P_137.bam',
        #'work/RNA_seq/align/mergedfinal/P_146.bam',
        #'work/RNA_seq/align/mergedfinal/P_254.bam',
        #'work/RNA_seq/align/mergedfinal/P_256.bam',
        #'work/RNA_seq/align/mergedfinal/P_257.bam',
        #'work/RNA_seq/align/mergedfinal/P_258.bam',
        #'work/RNA_seq/align/mergedfinal/P_265.bam',
        #'work/RNA_seq/align/mergedfinal/P_327.bam',
        #'work/RNA_seq/align/mergedfinal/P_266.bam',
        #'work/RNA_seq/align/mergedfinal/P_335.bam',
        #'work/RNA_seq/align/mergedfinal/P_348.bam',
        #'work/RNA_seq/align/mergedfinal/P_371.bam',
        #'work/RNA_seq/align/mergedfinal/P_389.bam',
        #'work/RNA_seq/align/mergedfinal/P_395.bam',
        #'work/RNA_seq/align/mergedfinal/P_400.bam',
        #expand('work/RNA_seq/align/mergedfinal/{id}.bam.bai', id = config["merged_patient_IDS"])


rule mergebam:
    input:
        T0 = lambda wildcards: expand(f'work/RNA_seq/align/T0/{config["SRR_merge_0"][wildcards.srr0]}.bam'),
        T1 = lambda wildcards: expand(f'work/RNA_seq/align/T1/{config["SRR_merge_1"][wildcards.srr1]}.bam'),
        T2 = lambda wildcards: expand(f'work/RNA_seq/align/T2/{config["SRR_merge_2"][wildcards.srr2]}.bam'), 
        T3 = lambda wildcards: expand(f'work/RNA_seq/align/T3/{config["SRR_merge_3"][wildcards.srr3]}.bam'), 
        T4 = lambda wildcards: expand(f'work/RNA_seq/align/T4/{config["SRR_merge_4"][wildcards.srr4]}.bam'),
    output:
        merged_temp = 'work/RNA_seq/align/mergedbam/{srr0}_{srr1}_{srr2}_{srr3}_{srr4}/{id}.bam'
    threads: 20
    resources:
        mem_mb = 6000,
        walltime = "8:00",
        disk_scratch = 30
    shell:
      '''
         software/samtools-1.15.1/samtools merge --threads {threads} -o {output.merged_temp} {input.T0} {input.T1} {input.T2} {input.T3} {input.T4}

      '''
rule merge3:
    input:
        T0 = lambda wildcards: expand(f'work/RNA_seq/align/T0/SRR12609780.bam'),
        T2 = lambda wildcards: expand(f'work/RNA_seq/align/T2/SRR12609783.bam'),
        T3 = lambda wildcards: expand(f'work/RNA_seq/align/T3/SRR12609784.bam'),
        T4 = lambda wildcards: expand(f'work/RNA_seq/align/T4/SRR12609787.bam')

    output:
        merged_temp = 'work/RNA_seq/align/mergedbam/SRR12609780_SRR12609783_SRR12609784_SRR12609787/P_3.bam'
    shell:
        '''
         software/samtools-1.15.1/samtools merge --threads {threads} -o {output.merged_temp} {input.T0} {input.T2} {input.T3} {input.T4}
         software/samtools-1.15.1/samtools index {output.merged_temp}
        '''

rule merge157:
    input:
        T1 = lambda wildcards: expand(f'work/RNA_seq/align/T1/SRR12609496.bam'),
        T2 = lambda wildcards: expand(f'work/RNA_seq/align/T2/SRR12609499.bam'),
        T3 = lambda wildcards: expand(f'work/RNA_seq/align/T3/SRR12609501.bam'),
        T4 = lambda wildcards: expand(f'work/RNA_seq/align/T4/SRR12609503.bam')

    output:
        merged_temp = 'work/RNA_seq/align/mergedbam/SRR12609496_SRR12609499_SRR12609501_SRR12609503/P_157.bam'
    shell:
        '''
         software/samtools-1.15.1/samtools merge --threads {threads} -o {output.merged_temp} {input.T1} {input.T2} {input.T3} {input.T4}
         software/samtools-1.15.1/samtools index {output.merged_temp}
        '''

rule movebam1:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609770_SRR12609772_SRR12609774_SRR12609775_SRR12609776/P_1.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_1.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam2:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609780_SRR12609783_SRR12609784_SRR12609787/P_3.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_3.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam3:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609800_SRR12609801_SRR12609803_SRR12609804_SRR12609805/P_13.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_13.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam4:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609808_SRR12609809_SRR12609812_SRR12609813_SRR12609814/P_29.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_29.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam5:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609817_SRR12609818_SRR12609820_SRR12609821_SRR12609822/P_36.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_36.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam6:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609394_SRR12609395_SRR12609397_SRR12609399_SRR12609401/P_75.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_75.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam7:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609405_SRR12609406_SRR12609409_SRR12609411_SRR12609412/P_76.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_76.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam8:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609414_SRR12609415_SRR12609418_SRR12609420_SRR12609422/P_83.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_83.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam9:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609425_SRR12609426_SRR12609429_SRR12609431_SRR12609433/P_84.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_84.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam10:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609434_SRR12609435_SRR12609438_SRR12609439_SRR12609440/P_85.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_85.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam11:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609855_SRR12609856_SRR12609859_SRR12609861_SRR12609863/P_87.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_87.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam12:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609449_SRR12609450_SRR12609453_SRR12609454_SRR12609457/P_92.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_92.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam13:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609489_SRR12609866_SRR12609869_SRR12609871_SRR12609872/P_99.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_99.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam14:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609475_SRR12609476_SRR12609480_SRR12609482_SRR12609484/P_105.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_105.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam15:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609875_SRR12609876_SRR12609880_SRR12609882_SRR12609883/P_117.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_117.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam16:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609885_SRR12609886_SRR12609890_SRR12609892_SRR12609894/P_124.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_124.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam17:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609897_SRR12609898_SRR12609901_SRR12609903_SRR12609905/P_127.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_127.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam18:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609910_SRR12609911_SRR12609914_SRR12609916_SRR12609918/P_136.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_136.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam19:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609922_SRR12609923_SRR12609926_SRR12609928_SRR12609930/P_137.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_137.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam20:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609934_SRR12609935_SRR12609938_SRR12609940_SRR12609942/P_146.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_146.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam21:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609496_SRR12609499_SRR12609501_SRR12609503/P_157.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_157.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam22:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609506_SRR12609507_SRR12609509_SRR12609511_SRR12609513/P_175.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_175.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam23:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609516_SRR12609517_SRR12609520_SRR12609522_SRR12609523/P_176.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_176.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam24:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609526_SRR12609527_SRR12609530_SRR12609531_SRR12609533/P_201.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_201.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam25:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609535_SRR12609536_SRR12609539_SRR12609541_SRR12609543/P_203.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_203.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam26:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609546_SRR12609547_SRR12609550_SRR12609551_SRR12609552/P_221.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_221.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam27:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609574_SRR12609575_SRR12609578_SRR12609580_SRR12609581/P_227.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_227.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam28:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609974_SRR12609975_SRR12609978_SRR12609980_SRR12609981/P_254.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_254.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam29:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609982_SRR12609983_SRR12609986_SRR12609988_SRR12609990/P_256.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_256.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam30:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609994_SRR12609995_SRR12609998_SRR12610000_SRR12610001/P_257.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_257.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam31:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610003_SRR12610004_SRR12610007_SRR12610009_SRR12610011/P_258.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_258.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam32:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610014_SRR12610015_SRR12610018_SRR12610020_SRR12610021/P_265.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_265.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam33:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610053_SRR12610054_SRR12610057_SRR12609587_SRR12609589/P_266.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_266.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam34:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609595_SRR12609596_SRR12609599_SRR12609601_SRR12609603/P_268.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_268.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam35:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609608_SRR12609609_SRR12609612_SRR12609614_SRR12609616/P_278.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_278.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam36:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609626_SRR12609627_SRR12609629_SRR12609634_SRR12609636/P_294.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_294.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam37:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609641_SRR12609642_SRR12609645_SRR12609647_SRR12609649/P_297.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_297.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam38:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609676_SRR12609677_SRR12609680_SRR12610058_SRR12610060/P_302.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_302.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam39:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610032_SRR12610033_SRR12610036_SRR12610037_SRR12610038/P_327.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_327.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam40:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610063_SRR12610064_SRR12610066_SRR12610068_SRR12610069/P_335.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_335.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam41:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610071_SRR12610072_SRR12610075_SRR12610077_SRR12610079/P_348.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_348.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam42:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609651_SRR12609652_SRR12609655_SRR12609656_SRR12609658/P_367.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_367.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam43:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610092_SRR12610093_SRR12610096_SRR12610098_SRR12610101/P_371.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_371.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam44:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610115_SRR12610116_SRR12610119_SRR12610121_SRR12610122/P_389.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_389.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam45:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610123_SRR12610124_SRR12610127_SRR12610129_SRR12610131/P_395.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_395.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam46:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12610148_SRR12610149_SRR12610152_SRR12609682_SRR12609684/P_400.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_400.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam47:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609687_SRR12609688_SRR12609691_SRR12609693_SRR12609694/P_432.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_432.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''
rule movebam48:
    input:
        bams = 'work/RNA_seq/align/mergedbam/SRR12609696_SRR12609697_SRR12609700_SRR12609702_SRR12609704/P_461.bam'
    output:
        merged = 'work/RNA_seq/align/mergedfinal/P_461.bam'
    threads: 1
    shell:
        '''
        mv {input.bams} work/RNA_seq/align/mergedfinal/
        '''

rule index_0:
    input:
        bam = 'work/RNA_seq/align/mergedfinal/{id}.bam'
    output:
        simple_index = 'work/RNA_seq/align/mergedfinal/{id}.bam.bai'
    threads: 4
    shell:
        '''
         software/samtools-1.15.1/samtools index {input.bam}
        '''