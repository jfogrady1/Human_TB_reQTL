from pathlib import Path,PurePath

wildcard_constraints:
    region = r'\w+',
    run = r'\w+'
    
config.setdefault('Run_name_0', 'DV_variant_calling')
#if 'binding_paths' in config:
 #   for path in config['binding_paths']:
  #      workflow.singularity_args += f' -B {path}'

cohort_samples = config['samples'] if 'glob' not in config['samples'] else glob_wildcards(config["bam_path_0"] + config["samples"]["glob"]).sample

#rule all:
 #   input:
  #      expand('work/RNA_seq/variant_call/T0/{region}.{preset}.vcf.gz', region = config["regions"], preset = config["GL_config"]),
   #     expand('{name}/{region}.{preset}.vcf.gz',name=config['Run_name_0'],region=config['regions'],preset=config['GL_config']),
    #    expand('work/RNA_seq/variant_call/T0/deepvariant/{sample}.{region}.g.vcf.gz',sample=cohort_samples, region = config["regions"], allow_missing=True)



#NOTE: may need to be updated if deepvariant changes its internal parameters.
def make_custom_example_arguments(model):
    match model:
        case 'RNA' | 'WES':
            return '--channels \'\' --split_skip_reads'
        case 'PACBIO':
            return '--add_hp_channel --alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "199" --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels "0.12"'
        case 'WGS':
            return '--channels "insert_size"'
        case _:
            return '--channels \'\' --split_skip_reads'

def get_regions(region):
    if region == 'all':
        return ''
    elif isinstance(config['regions'][region],str) and Path(config['regions'][region]).exists():
        return ' --regions ' + '"' + " ".join([l.strip() for l in open(config['regions'][region])]) + '"'
    else:
        return ' --regions ' + '"' + f'{" ".join(map(str,config["regions"][region]))}' + '"'

BAM_EXT = config.get('bam_index','.bai')

#error can be caused by corrupted file, check gzip -t -v
rule deepvariant_make_examples:
    input:
        reference = multiext(config['reference'],'','.fai'),
        bam = lambda wildcards: multiext(f'work/RNA_seq/align/T0/{config["samples"][wildcards.sample]}','','.bai')
    output:
        examples = temp('work/RNA_seq/variant_call/T0/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz'),
        gvcf = temp('work/RNA_seq/variant_call/T0/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz'),
        info = temp('work/RNA_seq/variant_call/T0/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz.example_info.json')
    params:
        examples = lambda wildcards, output: PurePath(output['examples']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model_args = make_custom_example_arguments(config['model']),
        regions = lambda wildcards: get_regions(wildcards.region) 
    threads: 20
    resources:
        mem_mb = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('make_examples',{}).get('walltime','4h'),
        scratch = "1G"
    priority: 50
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {input.reference[0]} \
        --include_med_dp \
        --reads {input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        --sample_name {wildcards.sample} \
        {params.regions} \
        {params.model_args} \
        --task {wildcards.N}
        '''

def get_checkpoint(model):
    match model:
        case 'PACBIO' | 'WGS':
            return f'/opt/models/{model.lower()}/model.ckpt'
        case 'hybrid':
            return '/opt/models/hybrid_pacbio_illumina/model.ckpt'
        case 'RNA':
            return config['checkpoint']

rule deepvariant_call_variants:
    input:
        expand('work/RNA_seq/variant_call/T0/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True)
    output:
        temp('work/RNA_seq/varaint_call/T0/deepvariant/intermediate_results_{sample}_{region}/call_variants_output.tfrecord.gz')
    params:
        examples = lambda wildcards,input: PurePath(input[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = get_checkpoint(config['model'])
    threads: config.get('resources',{}).get('call_variants',{}).get('threads',12)
    resources:
        mem_mb = config.get('resources',{}).get('call_variants',{}).get('mem_mb',1500),
        walltime = config.get('resources',{}).get('call_variants',{}).get('walltime','24h'),
        scratch = "1G"
    priority: 90
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/call_variants \
        --outfile {output} \
        --examples {params.examples} \
        --checkpoint {params.model}
        '''

rule deepvariant_postprocess:
    input:
        reference = multiext(config['reference'],'','.fai'),
        variants = rules.deepvariant_call_variants.output[0],
        gvcf = expand('work/RNA_seq/variant_call/T0/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True)
    output:
        vcf = multiext('work/RNA_seq/variant_call/T0/deepvariant/{sample}.{region}.vcf.gz','','.tbi'),
        gvcf = multiext('work/RNA_seq/variant_call/T0/deepvariant/{sample}.{region}.g.vcf.gz','','.tbi')
    params:
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz')
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('postprocess',{}).get('mem_mb',40000),
        walltime = config.get('resources',{}).get('postprocess',{}).get('walltime','4h'),
        scratch = "1G"
    priority: 100
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/postprocess_variants \
        --ref {input.reference[0]} \
        --infile {input.variants} \
        --outfile {output.vcf[0]} \
        --gvcf_outfile {output.gvcf[0]} \
        --nonvariant_site_tfrecord_path {params.gvcf} \
        --novcf_stats_report
        '''

def get_GL_config(preset):
    match preset:
        case 'DeepVariantWGS' | 'DeepVariant_unfiltered' | 'DeepVariantWES_MED_DP':
            return preset
        case _:
            return config['GL_config'][preset]

rule GLnexus_merge:
    input:
        gvcfs = expand('work/RNA_seq/variant_call/T0/deepvariant/{sample}.{region}.g.vcf.gz',sample=cohort_samples, region = config["regions"],allow_missing=True),
        tbi = expand('work/RNA_seq/variant_call/T0/deepvariant/{sample}.{region}.g.vcf.gz.tbi',sample=cohort_samples, region = config["regions"], allow_missing=True)
    output:
        multiext('work/RNA_seq/variant_call/T0/{region}.{preset}.vcf.gz','','.tbi')
    params:
        preset = lambda wildcards: get_GL_config(wildcards.preset),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024
    threads: config.get('resources',{}).get('merge',{}).get('threads',12),
    resources:
        mem_mb = config.get('resources',{}).get('merge',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('merge',{}).get('walltime','4h'),
        scratch = '50G'
    priority: 25
    container: config.get('containers',{}).get('GLnexus','docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.3')
    shell:
        '''
        /usr/local/bin/glnexus_cli \
        --dir $TMPDIR/GLnexus.DB \
        --config {params.preset} \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {input.gvcfs} |\
        bcftools view - |\
        bgzip -@ 4 -c > {output[0]} 

        tabix -p vcf {output[0]}
        '''
