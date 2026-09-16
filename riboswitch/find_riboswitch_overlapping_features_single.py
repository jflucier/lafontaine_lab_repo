import os
import argparse
import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

GB_EXT = ".gb"
# checking for genes with location +/- OVERLAP_MARGIN in nt
OVERLAP_MARGIN = 250


def find_overlapping_features(start, end, genbank, acc, o):
    if not os.path.exists(genbank):
        print(f"Genbank not found: {genbank}")
        sys.exit(0)
    else:
        print(f"Genbank found: {genbank}")
        # a = acc.split(".")[0]
        for record in SeqIO.parse(genbank, "genbank"):
            if record.name != acc:
                continue

            # print(f"{record.name} {acc}")
            overlapping_feat = find_features(record, start, end)
            output_features(genbank, acc, start, end, overlapping_feat, o)
            break
        # genome = SeqIO.read(genbank, "genbank")

        # all_feat.append({
        #     'start': start,
        #     'end': end,
        #     'target': a,
        #     'acc': acc,
        #     'features': overlapping_feat
        # })

        # l = len(all_feat)
        # print(f"Ouputting results for {acc}: {l} overlapping hits identified")
        # output_features(genbank,acc,start,end,overlapping_feat)


def output_features(genbank, acc, start, end, overlapping_feat, o):
    # Print the overlapping features
    with open(o, "a") as file1:
        if len(overlapping_feat) > 0:
            for feature in overlapping_feat:
                if 'gene' in feature["f"].qualifiers:
                    g = feature["f"].qualifiers['gene'][0]
                elif 'locus_tag' in feature["f"].qualifiers:
                    g = feature["f"].qualifiers['locus_tag'][0]
                else:
                    print(f"key gene or locus_tag not found in feature: {feature}")
                    sys.exit(0)
                if 'product' in feature["f"].qualifiers:
                    n = feature["f"].qualifiers['product'][0]
                elif 'db_xref' in feature["f"].qualifiers:
                    n = feature["f"].qualifiers['db_xref'][0]
                elif 'note' in feature["f"].qualifiers:
                    n = feature["f"].qualifiers['note'][0]
                else:
                    n = ""
                # xref = feature["f"].qualifiers['db_xref']
                # prod = feature["f"].qualifiers['product']
                loc = feature["t"] + ":" + str(feature["f"].location.start) + "-" + str(feature["f"].location.end)
                file1.write(f"{genbank}\t{acc}\t{start}\t{end}\t{g}\t{loc}\t{n}\n")
        else:
            file1.write(f"{genbank}\t{acc}\t{start}\t{end}\t\t\t\n")


def find_features(genome, start, end):
    # Find all features that overlap the range
    overlapping_features = []
    for feature in genome.features:
        # if feature.type == 'CDS' and feature.qualifiers['gene'][0] == 'ROBO2':
        # print(f"{feature.qualifiers.get('gene')}")
        if feature.type == 'gene':
            # if isinstance(feature.location, CompoundLocation):
            #ENSAPOG00000001575.1
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
    argParser.add_argument("-s", "--start", help="model tsv result file", required=True, type=int)
    argParser.add_argument("-e", "--end", help="output filepath", required=True, type=int)
    argParser.add_argument("-g", "--genbank", help="genbank dir", required=True)
    argParser.add_argument("-a", "--acc", help="accession", required=True)
    argParser.add_argument("-o", "--out", help="output file", required=True)

    args = argParser.parse_args()
    s = args.start
    e = args.end
    g = args.genbank
    a = args.acc
    o = args.out

    # if os.path.exists(o):
    #     os.remove(o)

    find_overlapping_features(s, e, g, a, o)
