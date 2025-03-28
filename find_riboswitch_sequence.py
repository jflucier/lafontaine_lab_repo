import os
import argparse
import sys
import gzip
import pyfastx
import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import complement, reverse_complement
from Bio.SeqFeature import FeatureLocation

FA_EXT = ".fa"
GB_EXT = ".gb"
# checking for genes with location +/- OVERLAP_MARGIN in nt
# TARGET_OFFSET = 500

def get_sequences(in_f, outpath, genome_base):
    print(f"Reading input file: {in_f}")

    res = pd.read_csv(
        in_f,
        sep='\t',
        header=None,
        names=['genome', 'chr', 'start', 'end', 'strand', 'gene']
    )
    print(f"Grouping by genome path")
    x = res.groupby('genome')

    counter = 1
    out = []
    for gb_path, row in x:
        print(f"{counter}/{x.ngroups}: retreving sequences in {gb_path} hits")
        print(f"Total riboswitch to analyze: {len(row.index)}")
        counter += 1
        if not os.path.exists(gb_path):
            print(f"Genbank not found: {gb_path}")
            sys.exit(0)

        # try constructing fasta genome path
        # from: /fast2/def-lafontai/ensembl_genomes/genbank/anas_platyrhynchos/Anas_platyrhynchos.ASM874695v1.111.primary_assembly.5.dat
        # to: /fast2/def-lafontai/ensembl_genomes/fasta/anas_platyrhynchos/dna/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa
        # or
        # /fast2/def-lafontai/ensembl_genomes/fasta/anas_zonorhyncha/dna/Anas_zonorhyncha.ASM222487v1.dna.nonchromosomal.fa
        sp = os.path.basename(os.path.dirname(gb_path))
        base_fa_path = os.path.join(genome_base, "fasta", sp, "dna")
        pattern = re.compile(".*.nonchromosomal.dat$")
        if pattern.match(gb_path):
            # nonchromosomal
            fa_pattern = re.compile(".*.nonchromosomal.fa$")
        else:
            # toplevel
            fa_pattern = re.compile(".*.toplevel.fa$")

        fa_path = ""
        for f in os.listdir(base_fa_path):
            if fa_pattern.match(f):
                fa_path = os.path.join(base_fa_path,f)
                break

        if not os.path.exists(fa_path):
            print(f"Fasta not found. Falling back to toplevel file")
            fa_pattern = re.compile(".*.toplevel.fa$")
            for f in os.listdir(base_fa_path):
                if fa_pattern.match(f):
                    fa_path = os.path.join(base_fa_path, f)
                    break

        if not os.path.exists(fa_path):
            print(f"Fasta not found!")
            sys.exit(0)

        # g_records = SeqIO.to_dict(SeqIO.parse(gb_path, "genbank"))
        print(f"Loading genome in memory: {fa_path}")
        fa = pyfastx.Fasta(fa_path)
        # for name, seq in pyfastx.Fasta(fa_path, build_index=True):
        #     print(name, seq)

        # g_records = {}
        # with gzip.open(fa_path, "rt") as handle:
        #     for record in SeqIO.parse(handle, "fasta"):
        #         # print(record.accession)
        #         g_records[record.name] = record


        for _, r in row.iterrows():
            a = str(r['chr'])

            # if a in g_records:
            #     print(f"key {a} found")

            start = r['start']
            end = r['end']
            strand = r['strand']

            if strand == 1:
                seq = fa.fetch(a, (start, end))
            else:
                seq = reverse_complement(fa.fetch(a, (end, start)))

            out.append({
                'fa': fa_path,
                'chr': a,
                'start': start,
                'end': end,
                'strand': strand,
                'seq': seq
            })

        l = len(out)
        print(f"Outputting seq appended to lines")
        df = pd.DataFrame(out)
        df.to_csv(outpath, sep="\t", index=False)


def output_fa(out, all_seqs):
    with open(out, "w") as output:
        c = 1
        for f in all_seqs:
            print(">" + f['acc'] + "_" + str(c), file=output)
            print(f['seq'], file=output)
            c = c + 1


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()


    # mandatory
    argParser.add_argument("-i", "--input", help="input tsv: genome_file<tab>chr<tab>start<tab>end<tab>strand<tab>loc<tab>loc_strand<tab>overlap_type<tab>gene<tab>xref<tab>prod", required=True)
    argParser.add_argument("-o", "--output", help="output filepath", required=True)
    argParser.add_argument("-g", "--genome", help="base ensembl genome path filepath", required=True, default="/fast2/def-lafontai/ensembl_genomes/fasta")

    args = argParser.parse_args()
    i = args.input
    o = args.output
    g = args.genome

    get_sequences(i, o, g)
