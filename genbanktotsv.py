import re
import gzip
import sys


def parse_genbank_for_genes(genbank_content):
    gene_data = []
    current_gene_name = None
    current_note_lines = []  # To handle multi-line notes
    in_gene_block = False

    lines = genbank_content.splitlines()
    for line in lines:
        stripped_line = line.strip()

        # Check for the start of a new 'gene' feature block
        # A 'gene' line followed by its qualifiers
        if stripped_line.startswith("gene ") and not line[
            0].isspace():  # Check for 'gene' at the start of feature block
            if in_gene_block and current_gene_name:  # If we were in a gene block and have data, save it
                gene_data.append((current_gene_name, " ".join(current_note_lines).strip()))

            in_gene_block = True
            current_gene_name = None
            current_note_lines = []

            # Attempt to extract gene name directly from the 'gene' line
            match_gene_on_line = re.search(r'/gene="([^"]+)"', line)
            if match_gene_on_line:
                current_gene_name = match_gene_on_line.group(1)

        # Process lines within a 'gene' feature block
        if in_gene_block:
            # Extract gene name if it wasn't found on the 'gene' line itself
            if "/gene=" in line and current_gene_name is None:
                match = re.search(r'/gene="([^"]+)"', line)
                if match:
                    current_gene_name = match.group(1)

            # Extract note, handling multi-line notes
            if "/note=" in line:
                # If a new note starts, clear previous note parts
                current_note_lines = []
                # Check if the note is a single line (ends with ")
                if line.rstrip().endswith('"'):
                    note_match = re.search(r'/note="([^"]+)"', line)
                    if note_match:
                        current_note_lines.append(note_match.group(1).strip())
                else:  # Multi-line note starts
                    # Extract content after /note="
                    current_note_lines.append(line.split('/note="', 1)[1].strip())
            elif current_note_lines and stripped_line:  # Continuation of a multi-line note
                # Check if the line is part of the note continuation or a new qualifier
                if not stripped_line.startswith('/') and stripped_line != '"':
                    current_note_lines.append(stripped_line.rstrip('"'))
                elif stripped_line.startswith('/'):  # A new qualifier starts, so the note has ended
                    if current_gene_name:
                        gene_data.append((current_gene_name, " ".join(current_note_lines).strip()))
                    in_gene_block = False  # End of this gene's qualifiers
                    current_gene_name = None
                    current_note_lines = []  # Clear for next potential gene

    # After loop, add the last gene's data if it was being processed
    if in_gene_block and current_gene_name:
        gene_data.append((current_gene_name, " ".join(current_note_lines).strip()))

    return gene_data


# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python genbanktotsv.py <input_gz_file> <output_tsv_file>")
    sys.exit(1)

# Get file names from command-line arguments
gz_file_name = sys.argv[1]
output_tsv_name = sys.argv[2]

# Read the content from the .gz file
try:
    with gzip.open(gz_file_name, 'rt', encoding='utf-8') as f:
        file_content = f.read()
except FileNotFoundError:
    print(f"Error: The input file '{gz_file_name}' was not found.")
    print("Please make sure the .gz file exists and the path is correct.")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    sys.exit(1)

# Parse the content
extracted_genes = parse_genbank_for_genes(file_content)

# Write the extracted data to a TSV file
with open(output_tsv_name, "w", encoding='utf-8') as outfile:
    outfile.write("Gene Name\tNote\n")
    for gene_name, note in extracted_genes:
        outfile.write(f"{gene_name}\t{note}\n")

print(f"Extraction complete. Data saved to '{output_tsv_name}'")