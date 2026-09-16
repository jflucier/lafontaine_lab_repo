from Bio import SeqIO
import gzip
import sys
import os


def parse_genbank_with_biopython(gz_file_path):
    gene_data = []
    try:
        # Open the gzipped file and iterate through records
        with gzip.open(gz_file_path, 'rt', encoding='utf-8') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "gene":
                        gene_name = feature.qualifiers.get('gene', [''])[0]

                        # Notes can be tricky, they might be in 'note' or 'product' or other qualifiers
                        # Biopython collects them as lists, so we join them.
                        note = " ".join(feature.qualifiers.get('note', []))

                        # Sometimes the "note" is actually the product description, or other info
                        if not note:  # If no explicit 'note' qualifier, check 'product'
                            note = " ".join(feature.qualifiers.get('product', []))

                        if gene_name:  # Only add if a gene name was found
                            gene_data.append((gene_name, note))
    except FileNotFoundError:
        print(f"Error: The input file '{gz_file_path}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while parsing the file with Biopython: {e}")
        sys.exit(1)

    return gene_data


# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python your_script_name.py <input_gz_file> <output_tsv_file>")
    sys.exit(1)

# Get file names from command-line arguments
gz_file_name = sys.argv[1]
output_tsv_name = sys.argv[2]

print(f"Attempting to open with Biopython: {os.path.abspath(gz_file_name)}")

extracted_genes = parse_genbank_with_biopython(gz_file_name)

# Write the extracted data to a TSV file
with open(output_tsv_name, "w", encoding='utf-8') as outfile:
    outfile.write("Gene Name\tNote\n")
    for gene_name, note in extracted_genes:
        # Replace newlines and tabs within the note for proper TSV format
        cleaned_note = note.replace('\n', ' ').replace('\t', ' ')
        outfile.write(f"{gene_name}\t{cleaned_note}\n")

print(f"Extraction complete. Data saved to '{output_tsv_name}' using Biopython.")
