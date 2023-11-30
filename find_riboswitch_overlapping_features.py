import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

GB_EXT = ".gb"
# checking for genes with location +/- OVERLAP_MARGIN in nt
OVERLAP_MARGIN = 250

def find_overlapping_features(in_f, genome_path, out):
    res = pd.read_csv(in_f, sep='\t', index_col="accession.1")
    x = res.groupby('accession.1')

    for acc, row in x:
        print(f"{acc}")

        for _, r in row.iterrows():
            a = r['target name']
            start = r['target from coord']
            end = r['target to coord']
            gb_path = os.path.join(genome_path, a + GB_EXT)

            if not os.path.exists(gb_path):
                print(f"Genbank not found: {gb_path}")
                sys.exit(0)
            else:
                genome = SeqIO.read(gb_path, "genbank")
                overlapping_feat = find_features(genome, start, end)
                output_features(out, acc, r, overlapping_feat)


def output_features(out, acc, r, overlapping_feat):
    start = r['target from coord']
    end = r['target to coord']
    target = r['target name']
    # Print the overlapping features
    with open(out, "a") as output:
        for feature in overlapping_feat:
            g = feature.qualifiers['gene']
            xref = feature.qualifiers['db_xref']
            loc = feature.location
            output.write(f"{acc}\t{target}\t{start}\t{end}\t{g}\t{loc}\t{xref}\n")


def find_features(genome, start, end):
    # Find all features that overlap the range
    overlapping_features = []
    for feature in genome.features:
        if feature.type == 'gene':
            if isinstance(feature.location, FeatureLocation):
                if feature.location.start - OVERLAP_MARGIN <= end and feature.location.end + OVERLAP_MARGIN >= start:
                    overlapping_features.append(feature)

    return overlapping_features


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

    if os.path.exists(o):
        os.remove(o)

    find_overlapping_features(i, g, o)
