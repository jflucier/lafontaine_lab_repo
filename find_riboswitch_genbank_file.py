import argparse
import os
import re
import multiprocessing
from functools import partial
from tqdm import tqdm

def find_genbank_file(gen1, acc, file_index):
    """Finds a GenBank file using a pre-built index."""
    pattern = re.compile(rf"^{gen1}\..*\.{acc}\.dat$", re.IGNORECASE)

    for filename, filepath in file_index.items():
        if pattern.match(filename):
            return filepath
    return None

def find_nonchromosomal_file(gen1, file_index):
    """Finds a non-chromosomal GenBank file using regex."""
    pattern = re.compile(rf"^{gen1}\..*\.nonchromosomal\.dat$", re.IGNORECASE)

    for filename, filepath in file_index.items():
        if pattern.match(filename):
            return filepath
    return None

def process_line(line, gb_path, file_index):
    """Processes a single line and returns the output string."""
    fields = line.strip().split('\t')
    if len(fields) >= 11:
        gen1 = fields[0]
        acc = fields[2]
        s = fields[9]
        e = fields[10]
        g_file = find_genbank_file(gen1, acc, file_index)
        if g_file:
            return f"{s}\t{e}\t{acc}\t{g_file}\n"

        g_file = find_nonchromosomal_file(gen1, file_index)
        if g_file:
            return f"{s}\t{e}\t{acc}\t{g_file}\n"
    return None

def build_file_index(gb_path):
    """Builds an index of GenBank files for faster lookup."""
    file_index = {}
    for root, _, files in os.walk(gb_path):
        for file in files:
            file_index[file] = os.path.join(root, file)
    return file_index

def process_tsv_parallel(gb_path, input_file, output_file, num_processes=4):
    """Processes a TSV file in parallel with a progress bar and specified number of processes."""
    file_index = build_file_index(gb_path)

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines = infile.readlines()
            total_lines = len(lines)
            print(f"Processing {total_lines} lines using {num_processes} processes.")

            with multiprocessing.Pool(processes=num_processes) as pool:
                results = list(tqdm(pool.imap(partial(process_line, gb_path=gb_path, file_index=file_index), lines), total=total_lines))

            for result in results:
                if result:
                    outfile.write(result)

    except FileNotFoundError as e:
        print(f"Error: File not found: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TSV files and find corresponding GenBank files.")
    parser.add_argument("--rf_model", help="The RF model name.")
    parser.add_argument("--gb_path", help="The path to the GenBank directory.")
    parser.add_argument("--input_file", help="The input TSV file path.")
    parser.add_argument("--output_file", help="The output TSV file path.")
    parser.add_argument("--num_processes", type=int, default=4, help="Number of processes to use (default: 4).")

    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    if os.path.exists(args.output_file):
        os.remove(args.output_file)

    process_tsv_parallel(args.gb_path, args.input_file, args.output_file, args.num_processes)