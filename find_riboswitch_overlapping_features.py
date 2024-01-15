import os
import argparse
import sys

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
        print(f"Finding overlapping features for {acc} hits")
        all_feat = []
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
                overlapping_feat = find_features(genome, start, end)
                all_feat.append({
                    'start': start,
                    'end': end,
                    'target': a,
                    'acc': acc,
                    'features': overlapping_feat
                })

        l = len(all_feat)
        print(f"Ouputting results for {acc}: {l} overlapping hits identified")
        output_features(out, all_feat)


def output_features(out, all_feat):
    with open(out, "a") as output:
        for f in all_feat:
            start = f['start']
            end = f['end']
            target = f['target']
            acc = f['acc']
            overlapping_feat = f['features']
            # Print the overlapping features
            if len(overlapping_feat) > 0:
                for feature in overlapping_feat:
                    if 'gene' in feature["f"].qualifiers:
                        g = feature["f"].qualifiers['gene']
                    elif 'locus_tag' in feature["f"].qualifiers:
                        g = feature["f"].qualifiers['locus_tag']
                    else:
                        print(f"key gene or locus_tag not found in feature: {feature}")
                        sys.exit(0)
                    xref = feature["f"].qualifiers['db_xref']
                    loc = feature["t"] + ":" + str(feature["f"].location.start) + "-" + str(feature["f"].location.end)
                    output.write(f"{acc}\t{target}\t{start}\t{end}\t{g}\t{loc}\t{xref}\n")
            else:
                output.write(f"{acc}\t{target}\t{start}\t{end}\t\t\t\n")


def find_features(genome, start, end):
    # Find all features that overlap the range
    overlapping_features = []
    for feature in genome.features:
        # if feature.type == 'CDS' and feature.qualifiers['gene'][0] == 'ROBO2':
        if feature.type == 'CDS':
            # if isinstance(feature.location, CompoundLocation):
            if feature.location.start - OVERLAP_MARGIN <= end and feature.location.end + OVERLAP_MARGIN >= start:
                t = overlap_type(feature, start, end)
                overlapping_features.append({
                    "f": feature,
                    "t": t
                })

    return overlapping_features

def overlap_type(feature, start, end):
    if feature.location.start - OVERLAP_MARGIN >= start < feature.location.start:
        return 'UTR5'
    elif feature.location.end + OVERLAP_MARGIN <= end > feature.location.end:
        return 'UTR3'
    else:
        # loop through exons
        for l in feature.location.parts:
            if l.start >= start and l.end <= end:
                return 'exon'

        return 'intron'

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
