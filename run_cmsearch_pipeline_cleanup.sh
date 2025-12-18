#!/bin/bash

basepath=""
fasta_path=""
model_path=""
clean_index="-1"
TEMP=$(getopt -o b:f:m:i: --long base:,fa:,model:,index: -n 'rbs_submit.slurm' -- "$@")

eval set -- "$TEMP"

while true ; do
    case "$1" in
        -b|--base) basepath="$2" ; shift 2 ;;
        -f|--fa) fasta_path="$2" ; shift 2 ;;
        -m|--model) model_path="$2" ; shift 2 ;;
        -i|--index) clean_index="$2" ; shift 2 ;;
        --) shift ; break ;;
        *) echo "Invalid Option" ; exit 1 ;;
    esac
done

if [ "$clean_index" == "-1" ]; then
    echo "Error: --index must be specified."
    exit 1
fi

if [[ ! "$clean_index" =~ ^[0-9]+$ ]]; then
    echo "Error: --index '$clean_index' is not a valid integer."
    exit 1
fi

if [ ! -d "$fasta_path" ]; then
    echo "Error: Fasta path '$fasta_path' does not exist."
    exit 1
fi

if [ ! -d "$basepath" ]; then
    echo "Error: Base path '$basepath' does not exist."
    exit 1
fi

echo "Running the riboswitch pipeline"
echo "BASE_PATH: ${basepath}"
echo "FA_PATH: ${fasta_path}"
echo "MODEL_PATH: ${model_path}"

cd ${basepath}

export rf_model_path=${model_path}
export rf_model_file=$(basename $rf_model_path)
export rf_model="${rf_model_file%.cm}"
export b=$(basename ${basepath})

echo "running ${rf_model}"

outpath=${basepath}/${rf_model}
mkdir -p ${outpath}

echo "running cmsearch"
COUNTER=1
total=$(find "${fasta_path}" -name "*.toplevel.fa.gz" -not -path "*/dna_index/*" | wc -l)
find "${fasta_path}" -name "*.toplevel.fa.gz" -not -path "*/dna_index/*" -print0 | while IFS= read -r -d '' g; do
  genome_tmp=$(basename $g)
  gf=${genome_tmp%.gz}

  if [ "$COUNTER" -gt "$CLEAN_INDEX" ]; then
    echo "+++++ cleaning $COUNTER/${total}: ${outpath}/${gf}.${rf_model}.tsv"
    rm ${outpath}/${gf}.${rf_model}.tsv
  else
    echo "+++++ $COUNTER/${total}: results already exists: ${outpath}/${gf}.${rf_model}.tsv"
  fi

  COUNTER=$[$COUNTER +1]
done

echo "rm cmsearch results"
rm -f ${outpath}/all.${rf_model}.tsv
