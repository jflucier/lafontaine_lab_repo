import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

GB_EXT = ".gb"
def find_overlapping_features(in_f,genome_path,out):
    res = pd.read_csv(in_f, sep='\t', index_col="accession.1")
    x = res.groupby('accession.1')

    for acc, row in x:
        print(f"{acc}")

        for a in row['target name']:
            if not os.path.exists(os.path.join(genome_path, a + GB_EXT)):
                ## die
            else:
                genome = SeqIO.read(os.path.join(genome_path, a + GB_EXT), "genbank")
                start = 100
                end = 200
                overlapping_feat = find_features(genome,start,end)
                output_features(out,acc,a,overlapping_feat)


def output_features(out,acc,target,overlapping_feat):
    # Print the overlapping features
    with open(out, "a") as output:
        for feature in overlapping_feat:
            print(f"{acc}\t{target}\t{feature}")

def find_features(genome,start,end):
    # Find all features that overlap the range
    overlapping_features = []
    for feature in genome.features:
        if isinstance(feature.location, FeatureLocation):
            if feature.location.start <= end and feature.location.end >= start:
                overlapping_features.append(feature)

    return overlapping_features

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-i", "--input", help="model tsv result file", required=True)
    argParser.add_argument("-o", "--output", help="output filepath", required=True)
    argParser.add_argument("-g", "--genbank", help="genbank dir", required=True)

    args = argParser.parse_args()
    # i = args.input
    # o = args.output
    # g = args.genbank
    i = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/all.found.taxid.tsv"
    o = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/all.found.taxid.overlapping.tsv"
    g = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/refgenomes_hits"

    find_overlapping_features(i,g,o)