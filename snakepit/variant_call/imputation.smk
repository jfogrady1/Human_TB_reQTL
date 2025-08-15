autosomes = [str(i) for i in range(1, 23)]
autosome = [str(i) for i in range(1, 23)]
#if 'binding_paths' in config:
 #   for path in config['binding_paths']:
  #      workflow.singularity_args += f' -B {path}'


#rule all:
    #input:
        #'work/DNA_seq/imputation/all.filtered.vcf.gz',
        #expand('work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz.tbi', chromosome=autosomes),
        #expand('work/DNA_seq/imputation/all.filtered.chr{chromosomes}.bcf', chromosomes=autosomes),
        #expand('work/DNA_seq/imputation/chr_{chromosome}_phased.bcf', chromosome=autosomes),
        #expand('work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz', chromosome=autosomes),
        #'work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz',
        #'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz',
        #'work/DNA_seq/imputation/final/DV_IMPUTED_R0.3_filtered.vcf.gz'

rule filter:
    input:
        raw = 'work/RNA_seq/variant_call/mergedcall/all.Unrevised.vcf.gz',
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'

    output:
        filtered =  multiext('work/DNA_seq/imputation/all.filtered.vcf.gz','','.tbi')

    shell:
        '''
        {input.bcftool_exc} +fill-tags {input.raw} -Ov | {input.bcftool_exc} filter -e 'QUAL < 20' -Ov | {input.bcftool_exc} view -q 0.03:minor  -Ov | {input.bcftool_exc} filter -e 'HWE < 0.00001' | {input.bcftool_exc} filter -i 'F_MISSING < 0.2' -Oz -o {output.filtered[0]} && software/htslib-1.15.1/tabix -p vcf {output.filtered[0]}
        '''


rule split:
    input:
        filtered = 'work/DNA_seq/imputation/all.filtered.vcf.gz',
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'

    output:
        split = multiext('work/DNA_seq/imputation/all.filtered.chr{chromosomes}.bcf','','.csi')

    shell:
        '''
        {input.bcftool_exc} view -r chr{wildcards.chromosomes} {input.filtered} -Ob -o {output.split[0]} && {input.bcftool_exc} index {output.split[0]}
        '''


rule convert:
     input:
        vcf = 'data/1kg/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chromosomes}.filtered.shapeit2-duohmm-phased.vcf.gz',
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'

     output:
        bcf = multiext('data/1kg/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chromosomes}.filtered.shapeit2-duohmm-phased.bcf','','.csi')

     shell:
        '''
        {input.bcftool_exc} view {input.vcf} -Ob -o {output.bcf[0]} && {input.bcftool_exc} index {output.bcf[0]}
        '''


rule phase:
    input:
        target = multiext('work/DNA_seq/imputation/all.filtered.chr{chromosome}.bcf','','.csi'),
        reference = multiext('data/1kg/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chromosome}.filtered.shapeit2-duohmm-phased.bcf','','.csi'),
        map = 'data/1kg/chr{chromosome}.b38.gmap.gz'

    output:
        phased =  'work/DNA_seq/imputation/chr_{chromosome}_phased.bcf'

    container: 'software/shapeit5/shapeit5.sif'

    threads: 8
    shell:
        '''
        phase_common_static --input {input.target[0]} --reference {input.reference[0]} --region chr{wildcards.chromosome} --map {input.map} --thread {threads} --output {output.phased}
        ''' 

rule back2vcf: 
    input:
        bcf = 'work/DNA_seq/imputation/chr_{chromosome}_phased.bcf',
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'

    output:
        vcf = 'work/DNA_seq/imputation/chr_{chromosome}_phased.vcf.gz'

    shell: 
        '''
        {input.bcftool_exc} view {input.bcf} -Oz -o {output.vcf}
        '''      


rule imputation:
    input:
        gt= 'work/DNA_seq/imputation/chr_{chromosome}_phased.vcf.gz',
        ref =  'data/1kg/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz'

    output:
        out = 'work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz'

    threads: 8

    shell:
        '''
        java -jar software/beagle5.4/beagle.05May22.33a.jar gt={input.gt} ref={input.ref} out=work/DNA_seq/imputation/DV_chr{wildcards.chromosome}_IMPUTED nthreads={threads} impute=true
        '''

rule index_impute:
    input:
        vcf='work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz'

    output:
        index='work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz.tbi'

    shell:
        '''
        software/htslib-1.15.1/tabix -p vcf {input.vcf}
        '''

rule merge:
    input:
        vcf = expand('work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz', chromosome = autosome),
        index = expand('work/DNA_seq/imputation/DV_chr{chromosome}_IMPUTED.vcf.gz.tbi', chromosome = autosome),
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'
    output:
        merged = 'work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz'

    shell:
        '''
        {input.bcftool_exc} concat {input.vcf} -Oz -o {output.merged}
        '''

rule index_merged:
    input:
        merged = 'work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz',
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'
    output: 
        index = 'work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf.gz.csi'

    shell:
        '''
        {input.bcftool_exc} index -f {input.merged}
        '''


rule filterimpute3:
    input:
        raw=multiext('work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf', '.gz', '.gz.csi'),
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'
    
    output:
        filtered= 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.3_filtered.vcf.gz'

    shell:
        '''
        {input.bcftool_exc} filter -i 'DR2>=0.3' {input.raw[0]} | {input.bcftool_exc} view -q 0.05:minor -Oz -o {output.filtered}
        software/htslib-1.15.1/tabix -p vcf {output.filtered}
        '''

rule filterimpute6:
    input:
        raw=multiext('work/DNA_seq/imputation/final/DV_raw_IMPUTED.vcf', '.gz', '.gz.csi'),
        bcftool_exc = 'software/bcftools-1.15.1/bcftools'

    output:
        filtered= 'work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz'

    shell:
        '''
        {input.bcftool_exc} filter -i 'DR2>=0.6' {input.raw[0]} | {input.bcftool_exc} view -q 0.05:minor -Oz -o {output.filtered}
        software/htslib-1.15.1/tabix -p vcf {output.filtered}
        '''

