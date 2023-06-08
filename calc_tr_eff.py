import os
import sys
import time
import json

import pandas as pd
import directory as DIR


runthreadN = 12


def calc_read_counts(seq_type):
    if seq_type == 'RF':
        in_star_dir = f'{DIR.STAR_BAM_RF_DIR}*.bam'
        out_star_file = f'{DIR.STAR_DIR}RF_read_counts.txt'
        in_bwa_dir = f'{DIR.BWA_BAM_RF_DIR}*.bam'
        out_bwa_file = f'{DIR.BWA_DIR}RF_read_counts.txt'
    else:
        in_star_dir = f'{DIR.STAR_BAM_TR_DIR}*.bam'
        out_star_file = f'{DIR.STAR_DIR}TR_read_counts.txt'
        in_bwa_dir = f'{DIR.BWA_BAM_TR_DIR}*.bam'
        out_bwa_file = f'{DIR.BWA_DIR}TR_read_counts.txt'
    gtf_file = DIR.SARS_COV2_GTF
    # -d -D
    os.system(f'featureCounts -a {gtf_file} -T {runthreadN} -M --fraction -o {out_star_file} {in_star_dir}')
    os.system(f'featureCounts -a {gtf_file} -T {runthreadN} -M --fraction -o {out_bwa_file} {in_bwa_dir}')

    return None


def calc_avg_read_counts(align):
    df_RF = pd.read_csv(f'{DIR.INITIAL_FILEPATH}{align}/RF_read_counts.txt', sep='\t', comment='#', index_col=0)
    df_RF = df_processing(df_RF)
    df_TR = pd.read_csv(f'{DIR.INITIAL_FILEPATH}{align}/TR_read_counts.txt', sep='\t', comment='#', index_col=0)
    df_TR = df_processing(df_TR)
    print(df_RF)
    print(df_TR)

    df_RF.to_csv(f'{DIR.INITIAL_FILEPATH}{align}/RF_avg_read_counts.txt', sep='\t')
    df_TR.to_csv(f'{DIR.INITIAL_FILEPATH}{align}/TR_avg_read_counts.txt', sep='\t')
    return None


def df_processing(df):
    df = df.rename(columns=lambda x: x.split('/')[-1])
    column_name = [col for col in df.columns if 'sorted' in col]
    df.drop(labels=column_name, axis=1, inplace=True)
    df = df.rename(columns=lambda x: x.split('.')[0])

    for hpi in DIR.HPI:
        hpi_column = [col for col in df.columns if f'hpi{hpi}_' in col]
        df[f'hpi{hpi}'] = 0
        for hpi_col in hpi_column:
            df[f'hpi{hpi}'] += df[hpi_col]
        df[f'hpi{hpi}'] /= len(hpi_column)
        df.drop(labels=hpi_column, axis=1, inplace=True)

    return df


def preprocess(seq_type):

    if seq_type == 'RF':
        iteration = 2
        bwa_bam_dir = DIR.BWA_BAM_RF_DIR
        star_bam_dir = DIR.STAR_BAM_RF_DIR
    else:
        iteration = 4
        bwa_bam_dir = DIR.BWA_BAM_TR_DIR
        star_bam_dir = DIR.STAR_BAM_TR_DIR
    for hpi in DIR.HPI:
        for i in range(iteration):
            # sam to bam (bwa)
            sam_to_bam(f'{DIR.BWA_DIR}hpi{hpi}/{seq_type}/{seq_type}{i + 1}.sam',
                       f'{bwa_bam_dir}hpi{hpi}_{seq_type}{i + 1}.bam')

            # copy bam (star)
            copy_bam(f'{DIR.STAR_DIR}hpi{hpi}/{seq_type}/{seq_type}{i+1}/alignments.bam',
                     f'{star_bam_dir}hpi{hpi}_{seq_type}{i+1}.bam')

            # bwa indexing and sort
            samtools_sort(f'{bwa_bam_dir}hpi{hpi}_{seq_type}{i + 1}.bam',
                          f'{bwa_bam_dir}hpi{hpi}_{seq_type}{i + 1}.sorted.bam')
            samtools_index(f'{bwa_bam_dir}hpi{hpi}_{seq_type}{i + 1}.sorted.bam')

            # star indexing and sort
            samtools_sort(f'{star_bam_dir}hpi{hpi}_{seq_type}{i + 1}.bam',
                          f'{star_bam_dir}hpi{hpi}_{seq_type}{i + 1}.sorted.bam')
            samtools_index(f'{star_bam_dir}hpi{hpi}_{seq_type}{i + 1}.sorted.bam')
    # end

    return None


def sam_to_bam(in_sam, out_bam):
    order = f'samtools view -Sb {in_sam} > {out_bam}'
    os.system(order)
    print(order)
    return None


def copy_bam(indir, outdir):
    os.system(f'cp {indir} {outdir}')
    return None


def samtools_sort(in_bam, out_bam):
    order = f'samtools sort -@ {runthreadN} -o {out_bam} {in_bam}'
    os.system(order)
    print(order)
    return None


def samtools_index(in_bam):
    order = f'samtools index -@ {runthreadN} {in_bam}'
    os.system(order)
    print(order)
    return None


MODE = sys.argv[1]
PARAM = sys.argv[2:]
if __name__ == '__main__':
    if MODE in locals().keys():
        locals()[MODE](*PARAM)
    else:
        sys.exit('error: cmd=%s' % MODE)
