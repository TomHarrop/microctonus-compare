#!/usr/bin/env python3

import pathlib2


#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))


###########
# GLOBALS #
###########

fasta_files = ['ma_FR_norm_k71_diplo1',
               'ma_IE_trim-decon_k71_diplo0',
               'ma_MA_norm_k71_diplo0',
               'mh_UNK_trim-decon_k41_diplo1']

# containers
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'


#########
# RULES #
#########

rule target:
    input:
        expand('output/busco/run_{fasta_file}/full_table_{fasta_file}.tsv',
               fasta_file=fasta_files)

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
