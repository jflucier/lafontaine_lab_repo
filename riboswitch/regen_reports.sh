#!/usr/bin/env bash
# regen_reports.sh -- regenerate the riboswitch report TSVs from the per-model sqlite DBs.
#
# Usage (from anywhere):
#   regen_reports.sh <run_dir> [gene_info.sqlite]
# e.g.
#   regen_reports.sh /path/ensembl_protists/20250715_riboswitch_protist_run1
#   (gene_info defaults to <run_dir>/../release-61/gene_info.sqlite)
#
# Expects <run_dir>/<rf_model>/<run>_<rf_model>.sqlite and writes, in <run_dir>/<rf_model>/ :
#   <rf_model>.ensembl.genome.riboswitch.report.tsv            (all hits, left join on gene_info)
#   <rf_model>.ensembl.genome.riboswitch.report.geneoverlap.tsv (only hits overlapping a gene)
set -uo pipefail

run_dir=$(realpath "${1:?usage: $0 <run_dir> [gene_info.sqlite]}")
gene_info=$(realpath "${2:-${run_dir}/../release-61/gene_info.sqlite}")
[[ -f "$gene_info" ]] || { echo "ERROR: gene_info not found: $gene_info" >&2; exit 1; }

# $1=sqlite file  $2=rf model (= table suffix)  $3=join type for gene_info (LEFT|INNER)  $4=output tsv
run_query() {
  local sqlitedb=$1 db_rf_model=$2 jointype=$3 out=$4
  local tmp="${out}.tmp"
  sqlite3 "$sqlitedb" '.separator "\t"' '.header on' "
ATTACH DATABASE '${gene_info}' AS db2;
select distinct
  rs.specie,
  rs.target_name chr,
  rs.query_name,
  rs.query_accession,
  rs.model_type,
  rs.model_from_coord,
  rs.model_to_coord,
  rs.target_from_coord,
  rs.target_to_coord,
  rs.strand_lbl,
  rs.trunc,
  rs.pass,
  rs.gc,
  rs.bias,
  rs.score,
  rs.evalue,
  rs.inc,
  rs.mdl_len,
  rs.seq_len,
  ann.start,
  ann.end,
  ann.strand,
  ann.gene ann_gene,
  g.id,
  g.info,
  g.gb_info,
  ann.loc,
  ann.loc_strand,
  ann.xref,
  ann.prod,
  ann.overlap_type,
  seq.seq
from
  rs_ensembl_genomes_${db_rf_model} rs
  left join rs_ensembl_genomes_annot_strandspec_${db_rf_model} ann
    on ann.specie=rs.specie
      AND ann.chr=rs.target_name
      AND ann.start=rs.target_from_coord
      AND ann.end=rs.target_to_coord
      AND ann.strand=rs.strand_lbl
  ${jointype} join db2.gene_info g on g.id=ann.gene
  left join rs_ensembl_genomes_annot_strandspec_seq_${db_rf_model} seq
    on seq.specie=rs.specie
      AND seq.chr=rs.target_name
      AND seq.start=rs.target_from_coord
      AND seq.end=rs.target_to_coord
      AND seq.strand=rs.strand_lbl
" > "$tmp" && mv "$tmp" "$out" || { rm -f "$tmp"; return 1; }
}

n_ok=0; n_fail=0
for d in "$run_dir"/RF*/; do
  rf_model=$(basename "$d")
  db_rf_model=$rf_model
  outpath=${d%/}
  sqlitedb=$(ls "$outpath"/*.sqlite 2>/dev/null | head -n1)
  if [[ -z "$sqlitedb" ]]; then
    echo "SKIP (no sqlite): $rf_model" >&2; n_fail=$((n_fail+1)); continue
  fi

  echo "generate report ${rf_model}.ensembl.genome.riboswitch.report.tsv"
  run_query "$sqlitedb" "$db_rf_model" LEFT  "${outpath}/${rf_model}.ensembl.genome.riboswitch.report.tsv" \
    || { echo "FAILED report: $rf_model" >&2; n_fail=$((n_fail+1)); continue; }

  echo "generate report ${rf_model}.ensembl.genome.riboswitch.report.geneoverlap.tsv"
  run_query "$sqlitedb" "$db_rf_model" INNER "${outpath}/${rf_model}.ensembl.genome.riboswitch.report.geneoverlap.tsv" \
    || { echo "FAILED geneoverlap: $rf_model" >&2; n_fail=$((n_fail+1)); continue; }

  n_ok=$((n_ok+1))
done
echo "done: ${n_ok} models ok, ${n_fail} failed/skipped"
[[ $n_fail -eq 0 ]]