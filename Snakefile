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

fasta_files = ['fopius_arisanus',
               'ma_FR_norm_k71_diplo1',
               'ma_IE_trim-decon_k71_diplo0',
               'ma_MA_norm_k71_diplo0',
               'mh_UNK_trim-decon_k41_diplo1']

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
mummer_container = 'shub://TomHarrop/singularity-containers:mummer_4.0.0beta2'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'

#########
# RULES #
#########

rule target:
    input:
        expand('output/busco/run_{fasta_file}/full_table_{fasta_file}.tsv',
               fasta_file=fasta_files),
        expand('output/stats/{fasta_file}.tsv',
               fasta_file=fasta_files),
        'output/plot_data/mummer_plot.jpg',
        'output/plot_data/mummer_test_data.Rds'

rule plot_alignments:
    input:
        plot_data = 'output/plot_data/mummer_combined.Rds'
    output:
        jpeg_file = 'output/plot_data/mummer_plot.jpg'
    log:
        'output/logs/plot_alignments.log'
    benchmark:
        'output/benchmarks/plot_alignments.tsv'
    threads:
        1
    singularity:
        r_container
    script:
        'src/plot_alignments.R'


rule make_distance_matrix:
    input:
        report_files = expand('output/mummer/{ref}-vs-{query}/out.report',
                              pairwise_combinations,
                              ref=fasta_files,
                              query=fasta_files)
    output:
        distance_matrix = 'output/plot_data/distance_matrix.Rds'
    log:
        'output/logs/make_distance_matrix.log'
    benchmark:
        'output/benchmarks/make_distance_matrix.tsv'
    threads:
        1
    singularity:
        r_container
    script:
        'src/distance_tree.R'

# plot data is a few hundred million points, make a small subset of data for testing the plot
rule subset_plot_data:  
    input:
        plot_data = 'output/plot_data/mummer_combined.Rds'
    output:
        plot_data = 'output/plot_data/mummer_test_data.Rds'
    log:
        'output/logs/subset_plot_data.log'
    benchmark:
        'output/benchmarks/subset_plot_data.tsv'
    threads:
        1
    singularity:
        r_container
    shell:
        'Rscript -e \"'
        'library(data.table) ; '
        'x <- readRDS(\'{input.plot_data}\') ; '
        'saveRDS(x[sample(.N, 10000)], \'{output.plot_data}\') \" '
        '&> {log}'

rule combine_mummer_plot_data:
    input:
        plot_data = expand('output/plot_data/mummer_{ref}-vs-{query}.Rds',
                           pairwise_combinations,
                           ref=fasta_files,
                           query=fasta_files)
    output:
        plot_data = 'output/plot_data/mummer_combined.Rds'
    log:
        'output/logs/combine_mummer_plot_data.log'
    benchmark:
        'output/benchmarks/combine_mummer_plot_data.tsv'
    threads:
        1
    singularity:
        r_container
    script:
        'src/combine_mummer_plot_data.R'

rule read_mummer_output:
    input:
        coords_file = 'output/mummer/{ref}-vs-{query}/out.1coords',
    output:
        plot_data = 'output/plot_data/mummer_{ref}-vs-{query}.Rds'
    log:
        'output/logs/read_mummer_output_{ref}-vs-{query}.log'
    benchmark:
        'output/benchmarks/read_mummer_output_{ref}-vs-{query}.tsv'
    threads:
        1
    singularity:
        r_container
    script:
        'src/read_mummer_output.R'

rule wga:
    input:
        ref = 'output/filtered_assemblies/{ref}_final-scaffolds.fa',
        query = 'output/filtered_assemblies/{query}_final-scaffolds.fa'
    output:
        'output/mummer/{ref}-vs-{query}/out.delta',
        'output/mummer/{ref}-vs-{query}/out.1coords',
        'output/mummer/{ref}-vs-{query}/out.report'
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
