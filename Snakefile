#!/usr/bin/env python3

import itertools
import pathlib2


#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))


def pairwise_combinations(x, y):
    '''
    Has to be output that snakemake likes. Each iteration of combo should
    yield something like: 
    (('ref', 'ma_FR_norm_k71_diplo1'), ('query', 'ma_MA_norm_k71_diplo0'))
    '''
    for combo in itertools.combinations(x, 2):
        yield(((combo[0]), (y[0][0], combo[1][1])))


###########
# GLOBALS #
###########

fasta_files = ['ma_FR_norm_k71_diplo1',
               'ma_IE_trim-decon_k71_diplo0',
               'ma_MA_norm_k71_diplo0',
               'mh_UNK_trim-decon_k41_diplo1']

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
mummer_container = 'shub://TomHarrop/singularity-containers:mummer_4.0.0beta2'

#########
# RULES #
#########

rule target:
    input:
        expand('output/busco/run_{fasta_file}/full_table_{fasta_file}.tsv',
               fasta_file=fasta_files),
        expand('output/stats/{fasta_file}.tsv',
               fasta_file=fasta_files),
        expand('output/mummer/{ref}-vs-{query}/out.delta',
               pairwise_combinations,
               ref=fasta_files,
               query=fasta_files)

rule wga:
    input:
        ref = 'output/filtered_assemblies/{ref}_final-scaffolds.fa',
        query = 'output/filtered_assemblies/{query}_final-scaffolds.fa'
    output:
        'output/mummer/{ref}-vs-{query}/out.delta'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                          'mummer_{ref}-vs-{query}.log'))
    benchmark:
        'output/benchmarks/mummer_{ref}-vs-{query}.tsv'
    params:
        wd = 'output/mummer/{ref}-vs-{query}',
        ref = lambda wildcards, input: resolve_path(input.ref),
        query = lambda wildcards, input: resolve_path(input.query)
    threads:
        20
    singularity:
        mummer_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'dnadiff {params.ref} {params.query} '
        '&> {log}'


rule filter_short_contigs:
    input:
        'data/microctonus_assemblies/{fasta_file}_final-scaffolds.fa'
    output:
        'output/filtered_assemblies/{fasta_file}_final-scaffolds.fa'
    benchmark:
        'output/benchmarks/filter_short_contigs_{fasta_file}.tsv'
    log:
        'output/logs/filter_short_contigs_{fasta_file}.log'
    threads:
        1
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'in={input} '
        'out={output} '
        'minlength=10000 '
        '2> {log}'


rule assembly_stats:
    input:
        'data/microctonus_assemblies/{fasta_file}_final-scaffolds.fa'
    output:
        'output/stats/{fasta_file}.tsv'
    benchmark:
        'output/benchmarks/stats_{fasta_file}.tsv'
    log:
        'output/logs/stats_{fasta_file}.log'
    threads:
        1
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input} '
        'minscaf=1000 '
        'format=3 '
        'threads={threads} '
        '> {output} '
        '2> {log}'

rule busco:
    input:
        fasta = 'data/microctonus_assemblies/{fasta_file}_final-scaffolds.fa',
        lineage = 'data/lineages/hymenoptera_odb9'
    output:
        'output/busco/run_{fasta_file}/full_table_{fasta_file}.tsv'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                          'busco_{fasta_file}.log'))
    benchmark:
        'output/benchmarks/busco_{fasta_file}.tsv'
    params:
        wd = 'output/busco',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage)
    threads:
        20
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--in {params.fasta} '
        '--out {wildcards.fasta_file} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species honeybee1 '
        '--mode genome '
        '&> {log}'
