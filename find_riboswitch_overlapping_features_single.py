import os
import argparse
import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

GB_EXT = ".gb"
# checking for genes with location +/- OVERLAP_MARGIN in nt
OVERLAP_MARGIN = 250

def find_overlapping_features(start, end, genbank, acc):
    if not os.path.exists(genbank):
        print(f"Genbank not found: {genbank}")
        sys.exit(0)
    else:
        print(f"Genbank found: {genbank}")
        # a = acc.split(".")[0]
        for record in SeqIO.parse(genbank, "genbank"):
            if record.name != acc:
                continue

            print(f"{record.name} {acc}")
            overlapping_feat = find_features(record, start, end)
            output_features(genbank, acc, start, end, overlapping_feat)
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


def output_features(genbank,acc,start,end,overlapping_feat):
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
            prod = feature["f"].qualifiers['product']
            loc = feature["t"] + ":" + str(feature["f"].location.start) + "-" + str(feature["f"].location.end)
            print(f"{genbank}\t{acc}\t{start}\t{end}\t{g}\t{loc}\t{xref}\t{prod}\n")
    else:
        print(f"{genbank}\t{acc}\t{start}\t{end}\t\t\t\n")



def find_features(genome, start, end):
    # Find all features that overlap the range
    overlapping_features = []
    for feature in genome.features:
        # if feature.type == 'CDS' and feature.qualifiers['gene'][0] == 'ROBO2':
        if feature.type == 'gene':
            # if isinstance(feature.location, CompoundLocation):
            if feature.location.start - OVERLAP_MARGIN <= int(end) and feature.location.end + OVERLAP_MARGIN >= int(start):
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
    argParser.add_argument("-s", "--start", help="model tsv result file", required=True)
    argParser.add_argument("-e", "--end", help="output filepath", required=True)
    argParser.add_argument("-g", "--genbank", help="genbank dir", required=True)
    argParser.add_argument("-a", "--acc", help="accession", required=True)

    args = argParser.parse_args()
    s = args.start
    e = args.end
    g = args.genbank
    a = args.acc

    # if os.path.exists(o):
    #     os.remove(o)

    find_overlapping_features(s, e, g, a)
