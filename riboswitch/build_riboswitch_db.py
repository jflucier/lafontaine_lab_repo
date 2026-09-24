#!/usr/bin/env python3
"""
build_riboswitch_db.py

Imports every <base>/ensembl_*/<run_dir>/<model>/<model>.ensembl.genome.riboswitch.report.tsv
into one SQLite table, prepending:
  ensembl_kingdom : genomes, fungi, plants, protists, metazoa  (from the ensembl_<kingdom> folder)
  run_id          : run folder name without the date and without a leading
                    riboswitch_ / rbs_ / ensemblgenomes_ prefix

The GTF-style "info" field is also split into info_gene_id, info_gene_source,
info_gene_biotype and info_other. Quotes inside fields are kept as-is.

Examples
  # full build
  python3 build_riboswitch_db.py /jbod2/def-lafontai/analysis/riboswitch_eukaryotes

  # quick test: 3 report files per (kingdom, run), max 200 rows per file -> *.test.sqlite
  python3 build_riboswitch_db.py BASE --per-run 3 --max-rows 200

  # quick test: only the first 10 files, only fungi
  python3 build_riboswitch_db.py BASE --limit 10 --kingdom fungi
"""
import argparse
import csv
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

TABLE = "riboswitch_report"
DEFAULT_BASE = "/jbod2/def-lafontai/analysis/riboswitch_eukaryotes"

# (name, type) in the order of the original report columns
REPORT_COLS = [
    ("specie", "TEXT"), ("chr", "TEXT"), ("query_name", "TEXT"),
    ("query_accession", "TEXT"), ("model_type", "TEXT"),
    ("model_from_coord", "INTEGER"), ("model_to_coord", "INTEGER"),
    ("target_from_coord", "INTEGER"), ("target_to_coord", "INTEGER"),
    ("strand_lbl", "INTEGER"), ("trunc", "TEXT"), ("pass", "INTEGER"),
    ("gc", "REAL"), ("bias", "REAL"), ("score", "REAL"), ("evalue", "REAL"),
    ("inc", "TEXT"), ("mdl_len", "INTEGER"), ("seq_len", "INTEGER"),
    ("start", "INTEGER"), ("end", "INTEGER"), ("strand", "INTEGER"),
    ("ann_gene", "TEXT"), ("id", "TEXT"), ("info", "TEXT"),
    ("gb_info", "TEXT"), ("loc", "TEXT"), ("loc_strand", "TEXT"),
    ("xref", "TEXT"), ("prod", "TEXT"), ("overlap_type", "TEXT"),
    ("seq", "TEXT"),
]
NCOL = len(REPORT_COLS)  # 32
INFO_IDX = [n for n, _ in REPORT_COLS].index("info")

# parsed out of the GTF-style "info" field: key "value"; key "value";
# keys not listed here are kept together in info_other
INFO_KEYS = ["gene_id", "gene_source", "gene_biotype"]
EXTRA_COLS = [f"info_{k}" for k in INFO_KEYS] + ["info_other"]
INFO_RE = re.compile(r'(\w+)\s+"([^"]*)"\s*;?')


def parse_info(info):
    if not info:
        return [None] * len(EXTRA_COLS)
    found, other = {}, []
    for k, v in INFO_RE.findall(info):
        if k in INFO_KEYS and k not in found:
            found[k] = v
        else:
            other.append(f'{k} "{v}"')
    return [found.get(k) for k in INFO_KEYS] + ["; ".join(other) or None]


def run_id_from_dir(run_dir):
    s = re.sub(r"^[0-9]{8}_", "", run_dir)
    s = re.sub(r"^(riboswitch|rbs)_", "", s)
    return re.sub(r"^ensemblgenomes_", "", s)


def select_files(base, args):
    files = sorted(base.glob("ensembl_*/*/*/*.ensembl.genome.riboswitch.report.tsv"))
    if args.kingdom:
        files = [f for f in files
                 if f.relative_to(base).parts[0].removeprefix("ensembl_") in args.kingdom]
    if args.run:
        files = [f for f in files if any(r in f.relative_to(base).parts[1] for r in args.run)]
    if args.per_run:
        seen, kept = defaultdict(int), []
        for f in files:
            key = f.relative_to(base).parts[:2]
            if seen[key] < args.per_run:
                kept.append(f)
                seen[key] += 1
        files = kept
    if args.limit:
        files = files[: args.limit]
    return files


def integrity_checks(con, expected_total, per_file_expected):
    """Print integrity checks; return number of failures."""
    fails = 0

    def report(label, ok, detail=""):
        nonlocal fails
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' - ' + detail) if detail else ''}")
        if not ok:
            fails += 1

    n = con.execute(f"SELECT COUNT(*) FROM {TABLE}").fetchone()[0]
    report("row count in DB == data rows read from files", n == expected_total,
           f"db={n} files={expected_total}")

    bad_groups = []
    for (k, r), exp in per_file_expected.items():
        got = con.execute(
            f"SELECT COUNT(*) FROM {TABLE} WHERE ensembl_kingdom=? AND run_id=?", (k, r)).fetchone()[0]
        if got != exp:
            bad_groups.append(f"{k}/{r}: db={got} files={exp}")
    report("per (kingdom, run) row counts match", not bad_groups, "; ".join(bad_groups))

    for name, typ in REPORT_COLS:
        if typ in ("INTEGER", "REAL"):
            bad = con.execute(
                f'SELECT COUNT(*) FROM {TABLE} WHERE "{name}" IS NOT NULL '
                f'AND typeof("{name}") NOT IN ("integer","real")').fetchone()[0]
            if bad:
                report(f"numeric column {name} holds only numbers", False, f"{bad} non-numeric values")
    numeric_ok = True  # failures for individual numeric columns are reported above
    bad = con.execute(f"SELECT COUNT(*) FROM {TABLE} WHERE strand_lbl NOT IN (1,-1) OR strand_lbl IS NULL").fetchone()[0]
    report("strand_lbl is always 1 or -1", bad == 0, f"{bad} rows" if bad else "")

    bad = con.execute(
        f"SELECT COUNT(*) FROM {TABLE} WHERE start IS NOT NULL AND "
        f"(start != target_from_coord OR end != target_to_coord OR strand != strand_lbl)").fetchone()[0]
    report("annotated rows: start/end/strand match the hit coordinates", bad == 0,
           f"{bad} rows" if bad else "")

    bad = con.execute(
        f"SELECT COUNT(*) FROM {TABLE} WHERE id IS NOT NULL AND info_gene_id IS NOT NULL AND id != info_gene_id").fetchone()[0]
    report("id == info_gene_id where both are set", bad == 0, f"{bad} rows" if bad else "")

    bad = con.execute(
        f"SELECT COUNT(*) FROM {TABLE} WHERE seq IS NOT NULL AND seq GLOB '*[^ACGTUNacgtun]*'").fetchone()[0]
    report("seq contains only nucleotide letters", bad == 0, f"{bad} rows" if bad else "")

    bad = con.execute(f"SELECT COUNT(*) FROM {TABLE} WHERE specie IS NULL OR chr IS NULL OR query_accession IS NULL").fetchone()[0]
    report("specie, chr, query_accession never empty", bad == 0, f"{bad} rows" if bad else "")
    return fails


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", nargs="?", default=DEFAULT_BASE, help="base directory (default: %(default)s)")
    ap.add_argument("--out", help="output sqlite path (default: <base>/riboswitch_eukaryotes.v1.sqlite, "
                                  "or riboswitch_eukaryotes.test.sqlite when a test flag is used)")
    ap.add_argument("--limit", type=int, help="TEST: only the first N report files")
    ap.add_argument("--per-run", type=int, metavar="N", help="TEST: at most N report files per (kingdom, run)")
    ap.add_argument("--max-rows", type=int, metavar="N", help="TEST: at most N data rows per file")
    ap.add_argument("--kingdom", action="append", help="only this kingdom (repeatable), e.g. --kingdom fungi")
    ap.add_argument("--run", action="append", help="only run dirs containing this text (repeatable)")
    ap.add_argument("--overwrite", action="store_true", help="replace the output db if it exists")
    ap.add_argument("--dry-run", action="store_true", help="list selected files and exit")
    args = ap.parse_args()

    base = Path(args.base)
    test_mode = any([args.limit, args.per_run, args.max_rows, args.kingdom, args.run])
    db = Path(args.out) if args.out else base / (
        "riboswitch_eukaryotes.test.sqlite" if test_mode else "riboswitch_eukaryotes.v1.sqlite")

    files = select_files(base, args)
    if not files:
        sys.exit("No report files found")
    if args.dry_run:
        for f in files:
            rel = f.relative_to(base)
            print(f"{rel.parts[0].removeprefix('ensembl_')}\t{run_id_from_dir(rel.parts[1])}\t{rel}")
        print(f"{len(files)} files")
        return

    if db.exists():
        if args.overwrite:
            db.unlink()
        else:
            sys.exit(f"ERROR: {db} already exists (use --overwrite or remove it)")

    con = sqlite3.connect(db)
    cols_sql = ",\n    ".join(
        ["ensembl_kingdom TEXT NOT NULL", "run_id TEXT NOT NULL"]
        + [f'"{n}" {t}' for n, t in REPORT_COLS]
        + [f'"{c}" TEXT' for c in EXTRA_COLS])
    con.execute(f"CREATE TABLE {TABLE} (\n    {cols_sql}\n)")
    insert_sql = f"INSERT INTO {TABLE} VALUES ({','.join('?' * (NCOL + 2 + len(EXTRA_COLS)))})"

    n_ok = n_skip = expected_total = 0
    per_file_expected = defaultdict(int)
    for f in files:
        rel = f.relative_to(base)
        kingdom = rel.parts[0].removeprefix("ensembl_")
        run_id = run_id_from_dir(rel.parts[1])
        rows, bad = [], None
        with open(f, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
            header = next(reader, None)
            if not header or header[0] != "specie":
                print(f"SKIP (unexpected header): {rel}", file=sys.stderr)
                n_skip += 1
                continue
            if len(header) != NCOL:
                print(f"WARN header has {len(header)} fields (expected {NCOL}): {rel}", file=sys.stderr)
            for lineno, row in enumerate(reader, start=2):
                if not row:
                    continue
                if len(row) > NCOL:
                    bad = f"{len(row)} fields (> {NCOL}) at line {lineno}"
                    break
                row += [""] * (NCOL - len(row))
                row = [v if v != "" else None for v in row]
                rows.append([kingdom, run_id] + row + parse_info(row[INFO_IDX]))
                if args.max_rows and len(rows) >= args.max_rows:
                    break
        if bad:
            print(f"SKIP ({bad}): {rel}", file=sys.stderr)
            n_skip += 1
            continue
        if not rows:
            print(f"SKIP (no data rows): {rel}", file=sys.stderr)
            n_skip += 1
            continue
        with con:
            con.executemany(insert_sql, rows)
        expected_total += len(rows)
        per_file_expected[(kingdom, run_id)] += len(rows)
        print(f"OK: {rel}  (kingdom={kingdom}, run_id={run_id}, rows={len(rows)})")
        n_ok += 1

    if n_ok == 0:
        con.close()
        db.unlink()
        sys.exit("No reports imported")

    print("\nIntegrity checks:")
    fails = integrity_checks(con, expected_total, per_file_expected)

    for name, cols in [("kingdom", "ensembl_kingdom"), ("run", "ensembl_kingdom, run_id"),
                       ("model", "query_accession"), ("specie", "specie"),
                       ("gene", "info_gene_id"), ("biotype", "info_gene_biotype"),
                       ("evalue", "evalue")]:
        con.execute(f"CREATE INDEX idx_{TABLE}_{name} ON {TABLE}({cols})")
    con.execute("ANALYZE")
    con.commit()

    print(f"\nDB: {db}")
    print(f"Imported {n_ok} reports ({expected_total} rows), skipped {n_skip}, integrity failures {fails}")
    print("kingdom\trun_id\tn\tannotated\twith_seq")
    for r in con.execute(f"SELECT ensembl_kingdom, run_id, COUNT(*), COUNT(start), COUNT(seq) "
                         f"FROM {TABLE} GROUP BY 1,2"):
        print("\t".join(map(str, r)))
    con.close()


if __name__ == "__main__":
    main()
