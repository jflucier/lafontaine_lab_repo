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
    genbank_dir_marker = os.path.join("genbank")  # Ensures cross-OS compatibility

    for gb_path, row in x:
        print(f"{counter}/{x.ngroups}: retreving sequences in {gb_path} hits")
        print(f"Total riboswitch to analyze: {len(row.index)}")
        counter += 1
        if not os.path.exists(gb_path):
            print(f"Genbank not found: {gb_path}")
            sys.exit(0)

        # try constructing fasta genome path
        # Use rsplit to split from the right, ensuring we get the correct last occurrence
        path_parts_split_by_genbank = gb_path.rsplit(genbank_dir_marker, 1)
        if len(path_parts_split_by_genbank) != 2:
            print(f"Error: Expected '{genbank_dir_marker}' not found or in unexpected place in path: {gb_path}")
            return None

        common_root_prefix = path_parts_split_by_genbank[0]  # e.g., "/fast2/def-lafontai/ensembl_protists/release-61/"
        specie_path = path_parts_split_by_genbank[1].lstrip('/')
        species_dir = os.path.dirname(specie_path)

        print(f"gb_path={gb_path}")
        print(f"common_root_prefix={common_root_prefix}")
        print(f"specie_path={specie_path}")
        print(f"species_dir={species_dir}")


        base_fa_path = os.path.join(common_root_prefix, "fasta", "fasta", species_dir, "dna")

        # from: /fast2/def-lafontai/ensembl_genomes/genbank/anas_platyrhynchos/Anas_platyrhynchos.ASM874695v1.111.primary_assembly.5.dat
        # to: /fast2/def-lafontai/ensembl_genomes/fasta/anas_platyrhynchos/dna/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa
        # or
        # /fast2/def-lafontai/ensembl_genomes/fasta/anas_zonorhyncha/dna/Anas_zonorhyncha.ASM222487v1.dna.nonchromosomal.fa

        # sp = os.path.basename(os.path.dirname(gb_path))
        # base_fa_path = os.path.join(genome_base, sp, "dna")
        print(f"base fa path: {base_fa_path}")
        pattern = re.compile(".*.dna.nonchromosomal.dat.gz$")
        if pattern.match(gb_path):
            # nonchromosomal
            fa_pattern = re.compile(".*.dna.nonchromosomal.fa.gz$")
        else:
            # toplevel
            fa_pattern = re.compile(".*.dna.toplevel.fa.gz$")

        fa_path = ""
        for root, dirs, files in os.walk(base_fa_path):
            for f in files:
                if fa_pattern.match(f):
                    fa_path = os.path.join(root, f)
                    break  # Break out of the inner loop (files)
            if fa_path:  # If fa_path was found, break out of the outer loop (os.walk)
                break

        if not os.path.exists(fa_path):
            print(f"Fasta not found. Falling back to toplevel file")
            fa_pattern = re.compile(".*.dna.toplevel.fa.gz$")
            for root, dirs, files in os.walk(base_fa_path):
                for f in files:
                    if fa_pattern.match(f):
                        fa_path = os.path.join(root, f)
                        break  # Break out of the inner loop (files)
                if fa_path:  # If fa_path was found, break out of the outer loop (os.walk)
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

            try:
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
            except Exception as e:
                print(f"######## Error fetching sequence {a} at {start}-{end}: {e} #######")
                continue

    if out:
        l = len(out)
        print(f"Outputting seq appended to lines. Nbr of line={l}")
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
