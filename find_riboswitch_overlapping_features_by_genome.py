import gzip
import os
import argparse
import sys
import re

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

GB_EXT = ".gb"
# checking for genes with location +/- OVERLAP_MARGIN in nt
OVERLAP_MARGIN = 250

def find_overlapping_features(in_f, out):
    print(f"Reading input file: {in_f}")
    res = pd.read_csv(in_f, sep='\t', index_col=3, header=None, names=['start','end','acc','genome'])
    print(f"Grouping by genome path")
    x = res.groupby('genome')

    counter = 1
    prev_gb_path = ""
    for gb_path, row in x:
        print(f"{counter}/{x.ngroups}: Finding overlapping features for {gb_path} hits")
        print(f"Total riboswitch to analyze: {len(row.index)}")
        counter += 1
        if not os.path.exists(gb_path):
            print(f"Genbank not found: {gb_path}")
            sys.exit(0)

        if prev_gb_path != gb_path:
            print(f"Loading genome in memory: {gb_path}")
            g_records_map = {}  # Clear previous records
            try:
                # Check if the file is gzipped and open accordingly
                if gb_path.endswith('.gz'):
                    with gzip.open(gb_path, "rt") as handle:  # 'rt' for read text mode
                        for record in SeqIO.parse(handle, "genbank"):
                            if record.id:
                                g_records_map[record.id] = record
                            if record.name and record.name != record.id:
                                g_records_map[record.name] = record
                else:
                    with open(gb_path, "r") as handle:
                        for record in SeqIO.parse(handle, "genbank"):
                            if record.id:
                                g_records_map[record.id] = record
                            if record.name and record.name != record.id:
                                g_records_map[record.name] = record
                prev_gb_path = gb_path
            except Exception as e:
                print(f"Error parsing Genbank file {gb_path}: {e}. Skipping this genome.")
                continue

        all_feat = []
        for _, r in row.iterrows():
            a = str(r['acc'])

            current_record = g_records_map.get(a)  # Safely get the record
            if current_record is None:
                print(f"Accession/Name '{a}' not found in loaded GenBank '{gb_path}'. Skipping.")
                continue

            start = int(r['start'])
            end = int(r['end'])

            if end < start:
                strand = -1
            else:
                strand = 1

            overlapping_feat = find_features(current_record, start, end, strand)
            all_feat.append({
                'start': start,
                'end': end,
                'strand': strand,
                'gb': gb_path,
                'acc': a,
                'features': overlapping_feat
            })

        l = len(all_feat)
        print(f"Outputting oerlapping features to {out}")
        output_features(out, all_feat)


def output_features(out, all_feat):
    with open(out, "a") as output:
        for f in all_feat:
            start = f['start']
            end = f['end']
            gb = f['gb']
            acc = f['acc']
            strand = str(f['strand'])
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
                    if 'db_xref' in feature["f"].qualifiers:
                        xref = feature["f"].qualifiers['db_xref']
                    else:
                        xref = ''
                    if 'product' in feature["f"].qualifiers:
                        prod = feature["f"].qualifiers['product']
                    else:
                        prod = ''
                    loc = feature["t"] + ":" + str(feature["f"].location.start) + "-" + str(feature["f"].location.end)
                    loc_strand = str(feature["f"].location.strand)
                    output.write(f"{gb}\t{acc}\t{start}\t{end}\t{strand}\t{g}\t{loc}\t{loc_strand}\t{xref}\t{prod}\n")


def find_features(genome, start, end, strand):
    # Find all features that overlap the range
    overlapping_features = []
    for feature in genome.features:
        if feature.type == 'CDS':
            if feature.location.strand == strand and feature.location.start - OVERLAP_MARGIN <= end and feature.location.end + OVERLAP_MARGIN >= start:
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
    # argParser.add_argument("-g", "--genbank", help="genbank dir", required=True)

    args = argParser.parse_args()
    i = args.input
    o = args.output
    if os.path.exists(o):
        os.remove(o)

    find_overlapping_features(i, o)
