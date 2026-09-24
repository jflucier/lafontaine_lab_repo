#!/bin/bash
# build_riboswitch_db.sh
# Usage: bash build_riboswitch_db.sh [base_dir]
#
# Imports every ensembl_*/<run_dir>/<model>/<model>.ensembl.genome.riboswitch.report.tsv
# into one SQLite table, prepending two columns:
#   ensembl_clade : genomes, fungi, plants, protists, metazoa
#   run_id        : run folder name without the date (and without riboswitch_/rbs_/ensemblgenomes_ prefix)
set -euo pipefail

BASE=${1:-/jbod2/def-lafontai/analysis/riboswitch_eukaryotes}
DB=${BASE}/riboswitch_eukaryotes.v1.sqlite
TABLE=riboswitch_report
NCOL=32   # number of columns in the original report
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

cd "$BASE"

if [[ -e "$DB" ]]; then
    echo "ERROR: $DB already exists, remove it first" >&2
    exit 1
fi

sqlite3 "$DB" <<EOF
CREATE TABLE $TABLE (
    ensembl_clade      TEXT NOT NULL,
    run_id             TEXT NOT NULL,
    specie             TEXT,
    chr                TEXT,
    query_name         TEXT,
    query_accession    TEXT,
    model_type         TEXT,
    model_from_coord   INTEGER,
    model_to_coord     INTEGER,
    target_from_coord  INTEGER,
    target_to_coord    INTEGER,
    strand_lbl         INTEGER,
    trunc              TEXT,
    pass               INTEGER,
    gc                 REAL,
    bias               REAL,
    score              REAL,
    evalue             REAL,
    inc                TEXT,
    mdl_len            INTEGER,
    seq_len            INTEGER,
    start              INTEGER,
    end                INTEGER,
    strand             INTEGER,
    ann_gene           TEXT,
    id                 TEXT,
    info               TEXT,
    gb_info            TEXT,
    loc                TEXT,
    loc_strand         TEXT,
    xref               TEXT,
    prod               TEXT,
    overlap_type       TEXT,
    seq                TEXT
);
EOF

n_ok=0
n_skip=0

for f in ensembl_*/*/*/*.ensembl.genome.riboswitch.report.tsv; do
    [[ -f "$f" ]] || continue

    kingdom=$(echo "$f" | cut -d/ -f1 | sed 's/^ensembl_//')   # genomes, fungi, ...
    run_dir=$(echo "$f" | cut -d/ -f2)
    run_id=$(echo "$run_dir" | sed -E 's/^[0-9]{8}_//; s/^(riboswitch|rbs)_//; s/^ensemblgenomes_//')

    if [[ $(wc -l < "$f") -lt 2 ]]; then
        echo "SKIP (no data rows): $f" >&2
        n_skip=$((n_skip+1))
        continue
    fi

    if [[ $(head -n1 "$f" | cut -f1) != "specie" ]]; then
        echo "SKIP (unexpected header): $f" >&2
        n_skip=$((n_skip+1))
        continue
    fi

    # sqlite3 .import parses double quotes even in tab mode; refuse rather than corrupt rows
    if grep -q '"' "$f"; then
        echo "SKIP (contains double quotes, import would be unsafe): $f" >&2
        n_skip=$((n_skip+1))
        continue
    fi

    # skip header, pad to NCOL fields, prepend kingdom + run_id
    if ! awk -F'\t' -v OFS='\t' -v clade="$kingdom" -v rid="$run_id" -v n="$NCOL" -v f="$f" '
        NR==1 { next }
        { sub(/\r$/, "") }
        NF > n { print "ERROR: " NF " fields (> " n ") in " f " line " NR > "/dev/stderr"; exit 1 }
        {
            for (i = NF+1; i <= n; i++) $i = ""
            print clade, rid, $0
        }
    ' "$f" > "$TMP/current.tsv"; then
        echo "SKIP (bad rows): $f" >&2
        n_skip=$((n_skip+1))
        continue
    fi

    sqlite3 "$DB" <<EOF
.mode tabs
.import '$TMP/current.tsv' $TABLE
EOF

    echo "OK: $f  (clade=$clade, run_id=$run_id)"
    n_ok=$((n_ok+1))
done

if [[ $n_ok -eq 0 ]]; then
    echo "No reports imported" >&2
    exit 1
fi

# empty strings -> NULL
cols=$(sqlite3 "$DB" "SELECT group_concat(name, ' ') FROM pragma_table_info('$TABLE');")
for c in $cols; do
    sqlite3 "$DB" "UPDATE $TABLE SET \"$c\" = NULL WHERE \"$c\" = '';"
done

sqlite3 "$DB" <<EOF
CREATE INDEX idx_${TABLE}_clade      ON $TABLE(ensembl_clade);
CREATE INDEX idx_${TABLE}_run        ON $TABLE(ensembl_clade, run_id);
CREATE INDEX idx_${TABLE}_model      ON $TABLE(query_accession);
CREATE INDEX idx_${TABLE}_specie     ON $TABLE(specie);
CREATE INDEX idx_${TABLE}_evalue     ON $TABLE(evalue);
ANALYZE;
EOF

echo "Imported $n_ok reports, skipped $n_skip"
sqlite3 -header "$DB" "SELECT ensembl_clade, run_id, COUNT(*) AS n, COUNT(start) AS annotated, COUNT(seq) AS with_seq FROM $TABLE GROUP BY 1,2;"
