import sys
import os
import time
import directory as DIR


runthreadN = 12     # you should check


def main_star(seq_type):
    if seq_type == 'RF':
        iteration = 2
    else:
        iteration = 4
    indexing_STAR(seq_type)
    for hpi in DIR.HPI:
        for i in range(iteration):
            print(time.ctime(), 'start alignment')
            filepath = f'{DIR.INITIAL_FILEPATH}star/hpi{hpi}/{seq_type}/{seq_type}{i + 1}/'
            os.system(f'mkdir -p {filepath}')
            infile = f'{DIR.INITIAL_FILEPATH}data/hpi{hpi}/{seq_type}/{seq_type}{i + 1}.reads.5.fastq'
            outdir = filepath
            alignment_STAR(infile, outdir, seq_type)
            print(time.ctime(), 'end alignment')
        # end preprocess for 1 hpi (RPF-seq files:2, mRNA-seq files:4)
    # end preprocess for all hpi

    return None


def indexing_STAR(seq_type):

    if seq_type == 'TR':
        star_index_order = f'STAR --runMode genomeGenerate --genomeDir {DIR.STAR_INDEX_DIR}TR/' \
                           f' --genomeFastaFiles {DIR.SARS_COV2_GENOME} --genomeSAindexNbases 6'
    else:
        star_index_order = f'STAR --runMode genomeGenerate --genomeDir {DIR.STAR_INDEX_DIR}RF/' \
                           f' --genomeFastaFiles {DIR.SARS_COV2_GENOME} --sjdbFileChrStartEnd {DIR.SARS_COV2_GENOME_ANNOTATION}' \
                           f' --sjdbOverhang 100 --genomeSAindexNbases 6'
    print(time.ctime(), star_index_order)
    os.system(star_index_order)
    return


def alignment_STAR(infile, outdir, seq_type):
    if seq_type == 'TR':
        index = 6
        star_index_dir = f'{DIR.STAR_INDEX_DIR}TR/'
    else:
        index = 1
        star_index_dir = f'{DIR.STAR_INDEX_DIR}RF/'
    out_bam_dir = outdir

    star_alignment_order = f'STAR --genomeDir {star_index_dir} --runThreadN {runthreadN}' \
                           f' --readFilesIn {infile} --limitIObufferSize 300000000 --limitOutSJcollapsed 10000000' \
                           f' --limitBAMsortRAM 10000000000 --outFileNamePrefix {out_bam_dir}' \
                           f' --outStd BAM_SortedByCoordinate --outReadsUnmapped Fastx' \
                           f' --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --outFilterType BySJout' \
                           f' --outFilterMismatchNoverReadLmax 0.04 --chimOutType WithinBAM HardClip' \
                           f' --alignEndsType EndToEnd --outSJfilterOverhangMin 12 12 12 12' \
                           f' --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1' \
                           f' --outSJfilterDistToOtherSJmin 0 0 0 0 --scoreGapNoncan 0 --scoreGapGCAG 0' \
                           f' --scoreGapATAC 0 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignSJoverhangMin 6' \
                           f' --alignSJDBoverhangMin {index} > {out_bam_dir}alignments.bam'

    print(time.ctime(), star_alignment_order)
    os.system(star_alignment_order)

    return None


def main_bwa(seq_type):
    if seq_type == 'RF':
        iteration = 2
    else:
        iteration = 4

    for hpi in DIR.HPI:
        for i in range(iteration):
            print(time.ctime(), 'start alignment')
            filepath = f'{DIR.INITIAL_FILEPATH}bwa/hpi{hpi}/{seq_type}/'
            os.system(f'mkdir -p {filepath}')
            infile = f'{DIR.INITIAL_FILEPATH}data/hpi{hpi}/{seq_type}/{seq_type}{i + 1}.reads.5.fastq'
            outfile = f'{filepath}{seq_type}{i+1}.sam'
            alignment_bwa(infile, outfile)
            print(time.ctime(), 'end alignment')
        # end preprocess for 1 hpi (RPF-seq files:2, mRNA-seq files:4)
    # end preprocess for all hpi
    return None


def indexing_bwa():
    bwa_index_order = f'bwa index {DIR.BWA_INDEX_DIR}sars_sequence.fa'
    print(time.ctime(), bwa_index_order)
    os.system(bwa_index_order)
    return None


def alignment_bwa(infile, outfile):
    bwa_alignment_order = f'bwa mem -t {runthreadN} {DIR.BWA_INDEX_DIR}sars_sequence.fa {infile} > {outfile}'
    print(time.ctime(), bwa_alignment_order)
    os.system(bwa_alignment_order)
    return None


# main_star('RF') main_star('TR')

# indexing_bwa()
# main_bwa('RF') main_bwa('TR')
MODE = sys.argv[1]
PARAM = sys.argv[2:]
if __name__ == '__main__':
    if MODE in locals().keys():
        locals()[MODE](*PARAM)
    else:
        sys.exit('error: cmd=%s' % MODE)