import os

def find_genbank_file(gb_path, gen1, acc, file_index):
    """Finds a GenBank file using a pre-built index."""
    key = f"{gen1}.*.{acc}.dat"
    if key in file_index:
        return file_index[key]
    return None

def find_nonchromosomal_file(gb_path, gen1, file_index):
    """Finds a non-chromosomal GenBank file using a pre-built index."""
    key = f"{gen1}.*.nonchromosomal.dat".lower()
    if key in file_index:
        return file_index[key]
    return None

def process_line(line, gb_path, file_index):
    """Processes a single line and returns the output string."""
    fields = line.strip().split('\t')
    if len(fields) >= 11:
        gen1 = fields[0]
        acc = fields[2]
        s = fields[9]
        e = fields[10]
        print(f"file_index: {file_index}")
        print(f"finding gb file using: {gb_path} - {gen1} - {acc}")
        g_file = find_genbank_file(gb_path, gen1, acc, file_index)
        if g_file:
            return f"{s}\t{e}\t{acc}\t{g_file}\n"

        g_file = find_nonchromosomal_file(gb_path, gen1, file_index)
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

def process_tsv_single_thread(rf_model, gb_path, input_file, output_file):
    """Processes a TSV file in single thread."""
    file_index = build_file_index(gb_path)  # Build the index once

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines = infile.readlines()
            total_lines = len(lines)
            print(f"Processing {total_lines} lines.")

            for line in lines:
                result = process_line(line, gb_path, file_index)
                if result:
                    outfile.write(result)

    except FileNotFoundError as e:
        print(f"Error: File not found: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    rf_model = "RF00174.DDC1"  # Replace with your rf_model
    gb_path = "/fast2/def-lafontai/ensembl_genomes/genbank"
    input_file = f"/fast2/def-lafontai/rf_run2/infernal/ensembl_genomes/{rf_model}.sp.test.tsv"
    output_file = f"/fast2/def-lafontai/rf_run2/infernal/ensembl_genomes/{rf_model}.sp.out2.tsv"

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if os.path.exists(output_file):
        os.remove(output_file)

    process_tsv_single_thread(rf_model, gb_path, input_file, output_file)