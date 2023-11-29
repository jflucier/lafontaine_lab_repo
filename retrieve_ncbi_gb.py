import os
import argparse
import pandas as pd
from Bio import Entrez
from Bio import SeqIO

RETMAX = 10**9
GB_EXT = ".gb"

def get_annotation(input_f, outdir, email, batchsize):

    def accessions_to_gb(accessions, db, batchsize, retmax):
        def batch(sequence, size):
            l = len(accessions)
            for start in range(0, l, size):
                yield sequence[start:min(start + size, l)]

        def extract_records(records_handle):
            genomes_it = SeqIO.parse(records_handle, "genbank")
            for record in genomes_it:
                yield record
            # genomes = [record for record in genomes_it]
            # count = 0
            # for record in genomes:
            #     with open(f"/home/jflucier/tmp/test{count}.gb", "w") as output_handle:
            #         SeqIO.write(record, output_handle, "gb")
            #         count = count + 1
            # buffer = []
            # for line in records_handle:
            #     if line.startswith(b'LOCUS') and buffer:
            #         # yield accession number and record
            #         yield buffer[0].split()[1], "".join(buffer)
            #         buffer = [line]
            #     else:
            #         buffer.append(line)
            # yield buffer[0].split()[1], "".join(buffer)

        def process_batch(accessions_batch):
            # get GI for query accessions
            query = " ".join(accessions_batch)
            query_handle = Entrez.esearch(db=db, term=query, retmax=retmax)
            gi_list = Entrez.read(query_handle)['IdList']

            # get GB files
            search_handle = Entrez.epost(db=db, id=",".join(gi_list))
            search_results = Entrez.read(search_handle)
            webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
            records_handle = Entrez.efetch(
                db=db,
                rettype="gbwithparts",
                retmax=batchsize,
                webenv=webenv,
                query_key=query_key,
                retmode="text"
            )
            yield from extract_records(records_handle)

        accession_batches = batch(accessions, batchsize)
        processed = 0
        tot = len(accessions)
        for acc_batch in accession_batches:
            processed = processed + batchsize
            if processed > tot:
                processed = tot
            print(f"Processing {processed} / {tot}")
            yield from process_batch(acc_batch)

    def write_record(dir, record):
        with open(os.path.join(dir, record.id + GB_EXT), "w") as output:
            SeqIO.write(record, output, "gb")

    res = pd.read_csv(input_f, sep='\t', index_col="accession.1")
    Entrez.email = email
    x = res.groupby('accession.1')

    for acc, row in x:
        print(f"{acc}")

        t_acc_lst = []
        for a in row['target name']:
            if not os.path.exists(os.path.join(outdir, a + GB_EXT))\
                    and a not in t_acc_lst:
                t_acc_lst.append(a)

        for record in accessions_to_gb(t_acc_lst, "nuccore", batchsize, RETMAX):
            write_record(outdir, record)

        print(f"Done {acc}")

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument("-i", "--input", help="model tsv result file", required=True)
    argParser.add_argument("-o", "--output", help="output dir", required=True)
    argParser.add_argument("-e", "--email", help="your entrez email", required=True)
    argParser.add_argument("-b", "--batchsize", nargs='?', help="download batch size (default 10)", const=10, type=int, default=10)

    args = argParser.parse_args()
    # f = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/all.found.taxid.tsv"
    # outdir = "/storage/Documents/service/biologie/lafontaine/20230920_riboswitch_eukaryotes/refgenomes_hits"
    # email = "jean-francois.lucier@usherbrooke.ca"
    # batch = 10
    get_annotation(
        args.input,
        args.output,
        args.email,
        args.batchsize
    )
