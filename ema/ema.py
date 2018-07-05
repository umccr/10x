#!/usr/bin/env python
import os
import sys
from glob import glob
from os.path import isfile, join, dirname, abspath
import shutil
import click
import subprocess

import yaml
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import verify_file, safe_mkdir, verify_dir
from ngs_utils import logger
from ngs_utils.logger import critical
from ngs_utils.utils import set_locale; set_locale()
from python_utils import hpc
from ngs_utils import snakemake_utils


def package_path():
    return dirname(abspath(__file__))


"""
ema.py /data/cephfs/punim0010/data/FASTQ/180312_A00130_0041_AHCLLMDMXX/Chromium_20180312/SI-GA-A11_*/*_R1_*.gz \
    -s NeverResponder10x_normal_EMA
    -g GRCh37
    -c
    -j 30
    -o ema_NeverResponder10x_normal_EMA
    
ema.py /data/cephfs/punim0010/data/FASTQ/180312_A00130_0041_AHCLLMDMXX/Chromium_20180312/SI-GA-A11_*/*_R1_*.gz \
    -s Colo829_100pc_10x
    -g GRCh37
    -c
    -j 30
"""


@click.command()
@click.argument('r1_fastq_paths', nargs=-1, type=click.Path(exists=True))
@click.option('-o', 'output_dir', type=click.Path(), help='Output directory (default is "ema_<sample_name>")')
@click.option('-j', '--jobs', 'jobs', default=1, help='Maximum number of cores to use at single time (works both for local '
              'and cluster runs)')
@click.option('-s', '--sample', 'sample_name', help='Sample name; required')
@click.option('--bins', help='Number of bins to split fastqs', default=500)
@click.option('-c', '--cluster-auto', 'cluster', is_flag=True, help='Submit jobs to cluster')
@click.option('-g', '--genome', help='Genome build (GRCh37 or hg38)', type=click.Choice(['GRCh37', 'hg38']), default='GRCh37')
@click.option('--unlock', is_flag=True, help='Propagaded to snakemake')
@click.option('--bc-whitelist', help='Whitelist of 10x barcodes')
@click.option('--bcbio-genomes', help='Path to bcbio-nextgen reference data (e.g. /bcbio/genomes; '
              'Hsapiens/{genome}/seq/{genome}.fa(.fai) and Hsapiens/{genome}/validation/giab-NA12878/truth_regions.bed are used)')
@click.option('--trim', '--trim-polyg', is_flag=True, help='Trim polyG in input reads')
def main(r1_fastq_paths, output_dir=None, jobs=None, sample_name=None, bins=None,
         cluster=False, genome=None, unlock=False, bc_whitelist=None, bcbio_genomes=None, trim_polyg=False):
    """
EMA wrapper.\n
r1_fastq_paths: paths to R1 fastq files for each sample\n
"""
    if not sample_name:
        critical('Please, specify sample name with -s')

    output_dir = output_dir or ('ema_' + sample_name)
    output_dir = safe_mkdir(abspath(output_dir))
    log_dir = safe_mkdir(join(output_dir, 'log'))
    logger.init(log_fpath_=join(log_dir, 'ema.log'), save_previous=True)

    #################################
    #### Making snakemake config ####
    #################################

    conf = dict()
    conf['package_path'] = package_path()
    conf['r1_fastq_paths'] = [verify_file(fp, is_critical=True) for fp in r1_fastq_paths]

    ###########################
    #### Setting ref paths ####
    ref_fa = None
    bwa_fa = None
    # Looking for reference fasta in bcbio_genomes folder?
    if bcbio_genomes:
        for subdir in [join(bcbio_genomes, f'Hsapiens/{genome}'),
                       join(bcbio_genomes, f'{genome}'),
                       join(bcbio_genomes)]:
            fp = join(subdir, 'seq', f'{genome}.fa')
            if isfile(fp):
                ref_fa = fp
            fp = join(subdir, 'bwa', f'{genome}.fa')
            if isfile(fp + '.bwt'):
                bwa_fa = fp
        if not ref_fa:
            critical(f'Not found seq/{genome}.fa in bcbio genomes directory {bcbio_genomes}')
        if not bwa_fa:
            critical(f'Not found BWA indices bwa/{genome}.fa* in bcbio genomes directory {bcbio_genomes}')

    # Reference files not provided, but we are on known host?
    if not ref_fa:
        ref_fa = hpc.get_ref_file(genome)
    if not bwa_fa:
        bwa_fa = hpc.get_ref_file(genome, ['bwa'], must_exist=False)
    if not bc_whitelist:
        bc_whitelist = hpc.get_loc().barcodes_10x

    # Done looking for refernece files, now some checks
    ref_fa = verify_file(ref_fa, is_critical=True, description='Reference fasta')
    verify_file(ref_fa + '.fai', is_critical=True, description='Reference fasta fai index')
    verify_file(bwa_fa + '.bwt', is_critical=True, description='BWA index')
    verify_file(bc_whitelist, is_critical=True)

    # Filling up conf
    conf['ref_fa'] = ref_fa
    conf['bwa_fa'] = bwa_fa
    conf['bc_whitelist'] = bc_whitelist

    #####################################
    #### Setting non-path parameters ####
    conf['sample'] = sample_name
    conf['bins'] = bins
    conf['trim_polyg'] = trim_polyg

    # Saving config to file
    conf_path = join(output_dir, 'conf.yml')
    with open(conf_path, 'w') as out:
        yaml.dump(conf, out)

    #########################
    #### Setting cluster ####
    #########################

    cluster_param = ''
    if cluster:
        log_dir = safe_mkdir(join(log_dir, 'cluster'))
        cluster_param = snakemake_utils.make_cluster_cmdl(log_dir)

    ###############################
    #### Building command line ####
    ###############################

    snakefile = join(package_path(), 'Snakefile')
    # import pdb; pdb.set_trace()
    cmd = (
        f'snakemake '
        f'--snakefile {snakefile} '
        f'--printshellcmds '
        f'--directory {output_dir} '
        f'-j {jobs} '
        f'--rerun-incomplete ' 
        f'{cluster_param} '
        f'--configfile {conf_path} '
        f'--restart-times 5 '
    )

    if unlock:
        print('* Unlocking previous run... *')
        run_simple(cmd + ' --unlock')
        print('* Now rerunning *')

    try:
        run_simple(cmd)
    except subprocess.CalledProcessError:
        logger.error('--------')
        logger.error(f'Error: snakemake returned a non-zero status. Working directory: {output_dir}')
        raise

    # logger.info('Generating snakemake benchmark report...')
    # run_simple(cmd + f' --report {log_dir}/snakemake_report.html')

    logger.error('--------')
    logger.info(f'Finished. Output directory: {output_dir}')

    logger.info('Cleaning up...', ending=' ')
    slurm_log_dir = safe_mkdir(join(log_dir, 'slurm'))
    for slurm_log in glob(f'{output_dir}/slurm-*.out'):
        shutil.move(slurm_log, slurm_log_dir)
    logger.info('Done.', print_date=False)


if __name__ == '__main__':
    main()
