import os, sys, time
import directory as DIR
# you should check alignment_to_human_genome because the code would not work if runthreadN is less than 10.


def main(seq_type):
    if seq_type == 'RF':
        iteration = 2
    else:
        iteration = 4
    adaptor_sequence = ['TGGAATTCTCGGGTGCCAAGG', 'GTTCAGAGTTCTACAGTCCGACGATC']
    for hpi in DIR.HPI:
        for i in range(iteration):
            filepath = f'{DIR.INITIAL_FILEPATH}data/hpi{hpi}/{seq_type}/{seq_type}{i+1}.'
            preprocess(filepath, adaptor_sequence)
        # end preprocess for 1 hpi (RPF-seq files:2, mRNA-seq files:4)
    # end preprocess for all hpi

    return None


def preprocess(filedir, adaptor_sequence):

    quality_trimming(filedir)
    adaptor_trimming(adaptor_sequence, filedir)
    remove_artifacts(filedir)
    remove_ncRNA_reads(filedir)
    remove_mycoplasma_reads(filedir)

    return None


def quality_trimming(filedir):
    infile = f'{filedir}reads.fastq'
    outfile = f'{filedir}reads1.fastq'
    order = f'{DIR.PYTHONPATH} {DIR.OPTDIR}hsptrim-1.2.6/hsptrim.py se -t sanger -d both -q 30 -l 17 -w 1 -e -5 {infile} {outfile}'
    print(time.ctime(), order)
    os.system(order)
    return None


def adaptor_trimming(adaptor_sequence, outdir):
    adaptor_3sequence = adaptor_sequence[0]
    adaptor_5sequence = adaptor_sequence[1]
    infile = f'{outdir}reads1.fastq'
    outfile = f'{outdir}reads2.fastq'
    order = f'cutadapt --quiet -O 10 -m 20 -a {adaptor_3sequence} -g {adaptor_5sequence} {infile} > {outfile}'
    print(time.ctime(), order)
    os.system(order)
    return None


def remove_artifacts(outdir):
    infile = f'{outdir}reads2.fastq'
    outfile = f'{outdir}reads3.fastq'
    order = f'{DIR.OPTDIR}fastx_toolkit-0.0.13/fastx_artifacts_filter -Q 33 -i {infile} -o {outfile}'
    print(time.ctime(), order)
    os.system(order)
    return None


def remove_ncRNA_reads(outdir):

    index_sequence = DIR.ncRNA_seq
    out_prefix = f'{outdir}sequence'
    infile = f'{outdir}reads3.fastq'
    outfile = f'{outdir}reads4.fastq'
    out_sam = f'{outdir}reads4.sam'
    out_log = f'{outdir}reads4.log'
    bowtie_index_order = f'{DIR.OPTDIR}bowtie2-2.1.0/bowtie2-build {index_sequence} {out_prefix}'
    bowtie_alignment_order = f'{DIR.OPTDIR}bowtie2-2.1.0/bowtie2 -k 1 --norc --very-sensitive -x {out_prefix}' \
                             f' --un {outfile} -S {out_sam} {infile} > {out_log}'
    print(time.ctime(), bowtie_index_order)
    os.system(bowtie_index_order)
    print(time.ctime(), bowtie_alignment_order)
    os.system(bowtie_alignment_order)
    return None


def remove_mycoplasma_reads(outdir):
    out_prefix = f'{outdir}sequence'
    infile = f'{outdir}reads4.fastq'
    outfile = f'{outdir}reads5.fastq'
    out_sam = f'{outdir}reads5.sam'
    out_log = f'{outdir}reads5.log'
    bowtie_index_order = f'{DIR.OPTDIR}bowtie2-2.1.0/bowtie2-build {DIR.MYCO_sequence} {out_prefix}'
    bowtie_alignment_order = f'{DIR.OPTDIR}bowtie2-2.1.0/bowtie2 -k 1 --very-sensitive -x {out_prefix} --un {outfile} -S {out_sam} {infile} > {out_log}'
    print(time.ctime(), bowtie_index_order)
    os.system(bowtie_index_order)
    print(time.ctime(), bowtie_alignment_order)
    os.system(bowtie_alignment_order)
    return None


# main('RF') main('TR)
MODE = sys.argv[1]
PARAM = sys.argv[2:]
if __name__ == '__main__':
    if MODE in locals().keys():
        locals()[MODE](*PARAM)
    else:
        sys.exit('error: cmd=%s' % MODE)

