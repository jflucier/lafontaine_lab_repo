import os
import argparse
import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

FA_EXT = ".fa"
GB_EXT = ".gb"
# checking for genes with location +/- OVERLAP_MARGIN in nt
TARGET_OFFSET = 500

def gen_multi_fasta(in_f, genome_path, out):
    res = pd.read_csv(in_f, sep='\t', index_col="accession.1")
    x = res.groupby('accession.1')


    for acc, row in x:
        print(f"Finding overlapping features for {acc} hits")
        all_seqs = []
        for _, r in row.iterrows():
            a = r['target name']
            start = r['target from coord']
            end = r['target to coord']
            gb_path = os.path.join(genome_path, a + GB_EXT)

            if not os.path.exists(gb_path):
                print(f"Genbank not found: {gb_path}")
                sys.exit(0)
            else:
                print(f"Genbank found: {gb_path}")
                genome = SeqIO.read(gb_path, "genbank")
                target_seq = get_sequence(genome, start, end)
                all_seqs.append({
                    'seq': target_seq,
                    'acc': a
                })

        fa_o = os.path.join(out, acc + FA_EXT)
        print(f"Ouputting multifasta for {acc} in {fa_o}")
        output_fa(fa_o, all_seqs)


def output_fa(out, all_seqs):
    with open(out, "w") as output:
        for f in all_seqs:
            print(">" + f['acc'], file=output)
            print(f['seq'], file=output)


def get_sequence(genome, start, end):
    # Find all features that overlap the range
    seq_start_idx = 0
    if start - TARGET_OFFSET < 0:
        seq_start_idx = 0
    else:
        seq_start_idx = start - TARGET_OFFSET

    seq_end_idx = len(genome.seq)
    if end + TARGET_OFFSET > seq_end_idx:
        seq_end_idx = len(genome.seq)
    else:
        seq_end_idx = end + TARGET_OFFSET

    return genome.seq[seq_start_idx:seq_end_idx]

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-i", "--input", help="model tsv result file", required=True)
    argParser.add_argument("-o", "--output", help="output filepath", required=True)
    argParser.add_argument("-g", "--genbank", help="genbank dir", required=True)

    args = argParser.parse_args()
    i = args.input
    o = args.output
    g = args.genbank
    # i = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/all.found.taxid.test.tsv"
    # o = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/all.found.taxid.overlapping.tsv"
    # g = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/refgenomes_hits"

    # if os.path.exists(o):
    #     os.remove(o)

    gen_multi_fasta(i, g, o)
