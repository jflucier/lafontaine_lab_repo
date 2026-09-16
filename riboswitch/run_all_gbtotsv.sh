#!/bin/bash

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <GENBANK_DIR> <TSV_DIR> <PYTHON_SCRIPT>"
    echo "Example: $0 /path/to/genbank/ /path/to/tsv/"
    exit 1
fi

# Assign command-line arguments to variables
GENBANK_DIR="$1"
TSV_DIR="$2"
PYTHON_SCRIPT="/home/def-lafontai/programs/lafontaine_lab_repo/genbanktotsv.py"

# Print the paths that will be used
echo "GenBank Directory: $GENBANK_DIR"
echo "TSV Output Directory: $TSV_DIR"
echo "Python Script Path: $PYTHON_SCRIPT"

# Create the output TSV directory if it doesn't exist
mkdir -p "$TSV_DIR"

# Get the total number of files to process for progress tracking
total_files=$(find "$GENBANK_DIR" -name "*.dat.gz" | wc -l)
echo "Total files to process: ${total_files}"

# Define the number of parallel jobs.
NUM_PARALLEL_JOBS=8 # <--- Adjust this number based on your system's CPU cores

echo "Starting parallel processing with ${NUM_PARALLEL_JOBS} jobs..."

# Find all .dat.gz files and pipe them to xargs for parallel processing
# We pass TSV_DIR and PYTHON_SCRIPT as arguments to bash -c
find "$GENBANK_DIR" -name "*.dat.gz" -print0 | xargs -0 -I {} -P "$NUM_PARALLEL_JOBS" bash -c '
  # Exit immediately if a command exits with a non-zero status.
  set -e

  # The arguments from xargs and the outer script are now available as $1, $2, $3...
  local_f="$1"         # This will be the current .gz file path from xargs
  local_tsv_dir="$2"   # This will be the TSV_DIR passed from the outer script
  local_python_script="$3" # This will be the PYTHON_SCRIPT path from the outer script

  b=$(basename "$local_f") # Extract just the filename (e.g., file.dat.gz)
  output_tsv_file="${local_tsv_dir}${b}.gb.tsv" # Construct the output TSV path

  # Get the directory and name of the Python script from its full path
  script_dir=$(dirname "$local_python_script")
  script_name=$(basename "$local_python_script")

  # Change working directory to where the Python script is located
  cd "$script_dir" || { echo "ERROR: Could not change to script directory $script_dir for file $b" >&2; exit 1; }

  echo "Processing: ${b}" # Indicate which file is currently being processed by this job

  # Run your Python script using its basename, relative to the *new* current directory
  # The arguments ($local_f and $output_tsv_file) remain absolute paths
  python "./$script_name" "$local_f" "$output_tsv_file" 2>&1

' _ {} "$TSV_DIR" "$PYTHON_SCRIPT" # Pass arguments: placeholder, file from xargs, TSV_DIR, PYTHON_SCRIPT

echo "All parallel jobs submitted. Waiting for completion..."
echo "All files processed."