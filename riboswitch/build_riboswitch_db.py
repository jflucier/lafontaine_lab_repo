#!/usr/bin/env python3
"""
build_riboswitch_db.py

Imports every <base>/ensembl_*/<run_dir>/<model>/<model>.ensembl.genome.riboswitch.report.tsv
into one SQLite table, prepending:
  ensembl_kingdom : name of the ensembl_<kingdom> folder (as-is, e.g. ensembl_fungi)
  run_id          : run folder name (as-is, e.g. 20250715_riboswitch_protist_run1)

Columns are matched BY HEADER NAME, not by position, so report files with a
different column set or order (e.g. no gb_info column) load correctly.
Columns missing from a file are stored as NULL.

The GTF-style "info" field is also split into info_gene_id, info_gene_source,
info_gene_biotype and info_other. Quotes inside fields are kept as-is.

Examples
  # full build
  python3 build_riboswitch_db.py /jbod2/def-lafontai/analysis/riboswitch_eukaryotes

  # quick test: 3 report files per (kingdom, run), max 200 rows per file -> *.test.sqlite
  python3 build_riboswitch_db.py BASE --per-run 3 --max-rows 200

  # quick test: only the first 10 files, only fungi
  python3 build_riboswitch_db.py BASE --limit 10 --kingdom fungi

  # list the distinct header layouts found, then exit
  python3 build_riboswitch_db.py BASE --scan-headers
"""
import argparse
import csv
import re
import sqlite3
import sys
from collections import defaultdict, Counter
from pathlib import Path

TABLE = "riboswitch_report"
DEFAULT_BASE = "/jbod2/def-lafontai/analysis/riboswitch_eukaryotes"

# canonical schema: (name, type) -- the DB column order
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
COL_NAMES = [n for n, _ in REPORT_COLS]
NCOL = len(REPORT_COLS)
REQUIRED_COLS = ["specie", "chr", "query_accession"]
ALLOWED_OVERLAP = ("intron", "exon", "UTR5", "UTR3")

# parsed out of the GTF-style "info" field: key "value"; key "value";
# keys not listed here are kept together in info_other
INFO_KEYS = ["gene_id", "gene_source", "gene_biotype"]
EXTRA_COLS = [f"info_{k}" for k in INFO_KEYS] + ["info_other"]
INFO_RE = re.compile(r'(\w+)\s+"([^"]*)"\s*;?')


def kingdom_from_dir(d):
    """ensembl_fungi -> fungi (works on Python < 3.9)"""
    return d[len("ensembl_"):] if d.startswith("ensembl_") else d


def run_id_from_dir(run_dir):
    """20250715_riboswitch_protist_run1 -> protist_run1 (not used by default)"""
    s = re.sub(r"^[0-9]{8}_", "", run_dir)
    s = re.sub(r"^(riboswitch|rbs)_", "", s)
    return re.sub(r"^ensemblgenomes_", "", s)


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


def select_files(base, args):
    files = sorted(base.glob("ensembl_*/*/*/*.ensembl.genome.riboswitch.report.tsv"))
    if args.kingdom:
        files = [f for f in files
                 if kingdom_from_dir(f.relative_to(base).parts[0]) in args.kingdom]
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


def read_header(f):
    """First line of a report as a list of column names ([] if the file is empty)."""
    if f.stat().st_size == 0:
        return []
    with open(f, newline="", encoding="utf-8") as fh:
        return next(csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE), [])


def scan_headers(files, base):
    """Group files by header layout and print a summary."""
    layouts = defaultdict(list)
    for f in files:
        layouts[tuple(read_header(f))].append(f.relative_to(base))
    print(f"{len(layouts)} distinct header layout(s) in {len(files)} files\n")
    for i, (hdr, fl) in enumerate(sorted(layouts.items(), key=lambda x: -len(x[1])), 1):
        missing = [n for n in COL_NAMES if n not in hdr]
        unknown = [h for h in hdr if h not in COL_NAMES]
        reordered = [h for h in hdr if h in COL_NAMES] != [n for n in COL_NAMES if n in hdr]
        print(f"layout {i}: {len(hdr)} columns, {len(fl)} files, e.g. {fl[0]}")
        print(f"   missing vs schema : {missing or '-'}")
        print(f"   unknown columns   : {unknown or '-'}")
        print(f"   order differs     : {reordered}")


def integrity_checks(con, expected_total, per_file_expected):
    """Print integrity checks; return number of failures."""
    fails = 0

    def report(label, ok, detail=""):
        nonlocal fails
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' - ' + detail) if detail else ''}")
        if not ok:
            fails += 1

    def count(where):
        return con.execute(f"SELECT COUNT(*) FROM {TABLE} WHERE {where}").fetchone()[0]

    n = count("1=1")
    report("row count in DB == data rows read from files", n == expected_total,
           f"db={n} files={expected_total}")

    bad_groups = []
    for (k, r), exp in per_file_expected.items():
        got = con.execute(
            f"SELECT COUNT(*) FROM {TABLE} WHERE ensembl_kingdom=? AND run_id=?", (k, r)).fetchone()[0]
        if got != exp:
            bad_groups.append(f"{k}/{r}: db={got} files={exp}")
    report("per (kingdom, run) row counts match", not bad_groups, "; ".join(bad_groups))

    bad_cols = []
    for name, typ in REPORT_COLS:
        if typ in ("INTEGER", "REAL"):
            bad = count(f"\"{name}\" IS NOT NULL AND typeof(\"{name}\") NOT IN ('integer','real')")
            if bad:
                bad_cols.append(f"{name}: {bad}")
    report("numeric columns hold only numbers", not bad_cols, "; ".join(bad_cols))

    bad = count("strand_lbl NOT IN (1,-1) OR strand_lbl IS NULL")
    report("strand_lbl is always 1 or -1", bad == 0, f"{bad} rows" if bad else "")

    bad = count("start IS NOT NULL AND (start != target_from_coord OR end != target_to_coord "
                "OR strand != strand_lbl)")
    report("annotated rows: start/end/strand match the hit coordinates", bad == 0,
           f"{bad} rows" if bad else "")

    bad = count("id IS NOT NULL AND info_gene_id IS NOT NULL AND id != info_gene_id")
    report("id == info_gene_id where both are set", bad == 0, f"{bad} rows" if bad else "")

    bad = count("seq IS NOT NULL AND seq GLOB '*[^ACGTUNMRWSYKVHDBacgtunmrwsykvhdb]*'")
    report("seq contains only nucleotide letters", bad == 0, f"{bad} rows" if bad else "")

    bad = count("start IS NOT NULL AND seq IS NULL")
    report("annotated rows always have a seq", bad == 0, f"{bad} rows" if bad else "")

    # catches column-shift problems: sequences / coordinates landing in a text column
    allowed = ",".join(f"'{v}'" for v in ALLOWED_OVERLAP)
    bad_vals = [r[0] for r in con.execute(
        f"SELECT DISTINCT substr(overlap_type,1,30) FROM {TABLE} "
        f"WHERE overlap_type IS NOT NULL AND overlap_type NOT IN ({allowed}) LIMIT 5")]
    report(f"overlap_type is NULL or one of {'/'.join(ALLOWED_OVERLAP)}", not bad_vals,
           f"e.g. {bad_vals}" if bad_vals else "")

    bad = count("loc_strand IS NOT NULL AND loc_strand NOT IN ('+','-','.','1','-1')")
    report("loc_strand looks like a strand symbol", bad == 0, f"{bad} rows" if bad else "")

    bad = count("specie IS NULL OR chr IS NULL OR query_accession IS NULL")
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
    ap.add_argument("--short-ids", action="store_true",
                    help="store kingdom as 'fungi' and run_id without date/prefix "
                         "(default: folder names as-is)")
    ap.add_argument("--overwrite", action="store_true", help="replace the output db if it exists")
    ap.add_argument("--dry-run", action="store_true", help="list selected files and exit")
    ap.add_argument("--scan-headers", action="store_true",
                    help="summarize the distinct header layouts of the selected files and exit")
    args = ap.parse_args()

    base = Path(args.base)
    test_mode = any([args.limit, args.per_run, args.max_rows, args.kingdom, args.run])
    db = Path(args.out) if args.out else base / (
        "riboswitch_eukaryotes.test.sqlite" if test_mode else "riboswitch_eukaryotes.v1.sqlite")

    files = select_files(base, args)
    if not files:
        sys.exit("No report files found")

    def ids(rel):
        if args.short_ids:
            return kingdom_from_dir(rel.parts[0]), run_id_from_dir(rel.parts[1])
        return rel.parts[0], rel.parts[1]

    if args.dry_run:
        for f in files:
            rel = f.relative_to(base)
            k, r = ids(rel)
            print(f"{k}\t{r}\t{rel}")
        print(f"{len(files)} files")
        return
    if args.scan_headers:
        scan_headers(files, base)
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
    insert_sql = f"INSERT INTO {TABLE} VALUES ({','.join('?' * (2 + NCOL + len(EXTRA_COLS)))})"

    n_ok = n_skip = n_empty = expected_total = 0
    per_file_expected = defaultdict(int)
    layout_counts = Counter()

    for f in files:
        rel = f.relative_to(base)
        kingdom, run_id = ids(rel)
        rows, bad = [], None

        if f.stat().st_size == 0:
            n_empty += 1
            continue

        with open(f, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
            header = next(reader, None)

            if not header:
                n_empty += 1
                continue

            if header[0] != "specie":
                print(f"SKIP (unexpected header): {rel}", file=sys.stderr)
                n_skip += 1
                continue

            unknown = [h for h in header if h not in COL_NAMES]
            missing_req = [c for c in REQUIRED_COLS if c not in header]
            dupes = [h for h, c in Counter(header).items() if c > 1]
            if unknown or missing_req or dupes:
                print(f"SKIP (unknown={unknown} missing_required={missing_req} "
                      f"duplicated={dupes}): {rel}", file=sys.stderr)
                n_skip += 1
                continue

            missing = [n for n in COL_NAMES if n not in header]
            if missing:
                print(f"NOTE {len(header)} columns, missing {missing} (stored as NULL): {rel}",
                      file=sys.stderr)
            layout_counts[tuple(header)] += 1

            # canonical column name -> index in THIS file (None if absent)
            idx = [header.index(n) if n in header else None for n in COL_NAMES]
            hlen = len(header)
            info_i = header.index("info") if "info" in header else None

            for lineno, row in enumerate(reader, start=2):
                if not row:
                    continue
                if len(row) > hlen:
                    bad = f"{len(row)} fields (> {hlen}) at line {lineno}"
                    break
                row += [""] * (hlen - len(row))
                vals = [row[i] if i is not None else "" for i in idx]
                vals = [v if v != "" else None for v in vals]
                info = row[info_i] if info_i is not None else None
                rows.append([kingdom, run_id] + vals + parse_info(info))
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

    print(f"\nHeader layouts imported: {len(layout_counts)}")
    for hdr, c in layout_counts.most_common():
        print(f"  {c} files, {len(hdr)} columns"
              + ("" if len(hdr) == NCOL else f", missing {[n for n in COL_NAMES if n not in hdr]}"))

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
    print(f"Imported {n_ok} reports ({expected_total} rows), empty/no-hit {n_empty}, "
          f"skipped errors {n_skip}, integrity failures {fails}")
    print("kingdom\trun_id\tn\tannotated\twith_seq")
    for r in con.execute(f"SELECT ensembl_kingdom, run_id, COUNT(*), COUNT(start), COUNT(seq) "
                         f"FROM {TABLE} GROUP BY 1,2"):
        print("\t".join(map(str, r)))
    con.close()


if __name__ == "__main__":
    main()
