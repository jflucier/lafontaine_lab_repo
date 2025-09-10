#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 30G
#SBATCH -A def-lafontai
#SBATCH --time 72:00:00
#SBATCH --output=rbs_%A.log
#SBATCH -J rbs

# Parsing command line options
ENV=""
basepath=""
fasta_path=""
gb_path=""
gene_info=""
model_path=""
TEMP=$(getopt -o e:b:f:g:i:m: --long env:,base:,fa:,gb:,gene:,model: -n 'rbs_submit.slurm' -- "$@")

eval set -- "$TEMP"

while true ; do
    case "$1" in
        -e|--env) ENV="$2" ; shift 2 ;;
        -b|--base) basepath="$2" ; shift 2 ;;
        -f|--fa) fasta_path="$2" ; shift 2 ;;
        -g|--gb) gb_path="$2" ; shift 2 ;;
        -i|--gene) gene_info="$2" ; shift 2 ;;
        -m|--model) model_path="$2" ; shift 2 ;;
        --) shift ; break ;;
        *) echo "Invalid Option" ; exit 1 ;;
    esac
done

if [[ "${ENV}" == "mp2" ]]; then
  echo "activating mp2 env"
  source /net/nfs-ip34/home/def-lafontai/programs/lafontaine_lab_repo/venv_mp2/bin/activate
  SCRIPT_PATH="/net/nfs-ip34/home/def-lafontai/programs/lafontaine_lab_repo/run_cmsearch_pipeline.sh"
else
  echo "activating ip34 env"
  source /home/def-lafontai/programs/lafontaine_lab_repo/venv/bin/activate
  SCRIPT_PATH="/home/def-lafontai/programs/lafontaine_lab_repo/run_cmsearch_pipeline.sh"
fi

sh ${SCRIPT_PATH} \
--base "${basepath}" \
--fa "${fasta_path}" \
--gb "${gb_path}" \
--gene "${gene_info}" \
--model "${model_path}"

echo "done!"
