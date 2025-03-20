import os
import subprocess

def process_tsv(rf_model, gb_path, input_file, output_file):
    """Processes a TSV file, finds corresponding GenBank files, and writes to an output file."""

    output_temp = f"infernal/ensembl_genomes/{rf_model}.sp.out.temp.tsv" #use temp file to avoid read/write conflicts

    try:
        with open(input_file, 'r') as infile, open(output_temp, 'w') as outfile:
            lines = infile.readlines()
            total_lines = len(lines)

            for counter, line in enumerate(lines, 1):
                print(f"running {counter} / {total_lines}")
                fields = line.strip().split('\t')
                if len(fields) >= 11:  # Ensure enough fields exist
                    gen1 = fields[0]
                    acc = fields[2]
                    s = fields[9]
                    e = fields[10]

                    g_file = find_genbank_file(gb_path, gen1, acc)

                    if g_file:
                        outfile.write(f"{s}\t{e}\t{acc}\t{g_file}\n")
                    else:
                        g_file = find_nonchromosomal_file(gb_path, gen1)
                        if g_file:
                            outfile.write(f"{s}\t{e}\t{acc}\t{g_file}\n")
                else:
                    print(f"Warning: Insufficient fields in line {counter}: {line.strip()}")

        os.rename(output_temp, output_file) #rename temp to final output.

    except FileNotFoundError as e:
        print(f"Error: File not found: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

def find_genbank_file(gb_path, gen1, acc):
    """Finds a GenBank file based on genome and accession."""
    for root, _, files in os.walk(gb_path):
        for file in files:
            if file == f"{gen1}.*.{acc}.dat":
                return os.path.join(root, file)
    return None

def find_nonchromosomal_file(gb_path, gen1):
    """Finds a non-chromosomal GenBank file based on genome."""
    for root, _, files in os.walk(gb_path):
        for file in files:
            if file.lower() == f"{gen1}.*.nonchromosomal.dat":
                return os.path.join(root, file)
    return None

if __name__ == "__main__":
    rf_model = "RF00174.DDC1"  # Replace with your rf_model
    gb_path = "/fast2/def-lafontai/ensembl_genomes/"
    input_file = f"/fast2/def-lafontai/rf_run2/infernal/ensembl_genomes/{rf_model}.sp.tsv"
    output_file = f"/fast2/def-lafontai/rf_run2/infernal/ensembl_genomes/{rf_model}.sp.out2.tsv"

    os.makedirs(os.path.dirname(output_file), exist_ok=True) #create output directory if it does not exist.
    if os.path.exists(output_file):
        os.remove(output_file)

    process_tsv(rf_model, gb_path, input_file, output_file)