#!/bin/bash
#source /home/def-lafontai/programs/lafontaine_lab_repo/venv/bin/activate
basepath=""
fasta_path=""
gb_path=""
gene_info=""
model_path=""
cmsearch=""
repo=""
tmp=""
TEMP=$(getopt -o b:f:g:i:m:c:s:t: --long base:,fa:,gb:,gene:,model:,cm:,script:,temp: -n 'rbs_submit.slurm' -- "$@")

eval set -- "$TEMP"

while true ; do
    case "$1" in
        -b|--base) basepath="$2" ; shift 2 ;;
        -f|--fa) fasta_path="$2" ; shift 2 ;;
        -g|--gb) gb_path="$2" ; shift 2 ;;
        -i|--gene) gene_info="$2" ; shift 2 ;;
        -m|--model) model_path="$2" ; shift 2 ;;
        -c|--cm) cm="$2" ; shift 2 ;;
        -s|--script) repo="$2" ; shift 2 ;;
        -t|--temp) tmp="$2" ; shift 2 ;;
        --) shift ; break ;;
        *) echo "Invalid Option" ; exit 1 ;;
    esac
done

echo "Running the riboswitch pipeline"
echo "BASE_PATH: ${basepath}"
echo "FA_PATH: ${fasta_path}"
echo "GB_PATH: ${gb_path}"
echo "GENE_INFO: ${gene_info}"
echo "MODEL_PATH: ${model_path}"
echo "CMSEARCH: ${cm}"
echo "REPO: ${repo}"
echo "TEMP: ${tmp}"

#basepath=/fast2/def-lafontai/ensembl_metazoa/20250902_riboswitch_metazoa_run2
#fasta_path=/fast2/def-lafontai/ensembl_metazoa/release-61/fasta/fasta/
#gb_path=/fast2/def-lafontai/ensembl_metazoa/release-61/genbank/genbank/
#gene_info=/fast2/def-lafontai/ensembl_metazoa/release-61/gene_info.sqlite
#model_path=/fast2/def-lafontai/ensembl_metazoa/20250902_riboswitch_metazoa_run2/updated_models

cd ${basepath}

export rf_model_path=${model_path}
export rf_model_file=$(basename $rf_model_path)
export rf_model="${rf_model_file%.cm}"
export db_rf_model="${rf_model//./_}"
export db_rf_model="${db_rf_model//-/_}"
export b=$(basename ${basepath})
export sqlitedb="${b}_${rf_model}.sqlite"

echo "running ${rf_model}"

outpath=${basepath}/${rf_model}
mkdir -p ${outpath}

echo "running cmsearch"
COUNTER=1
#rm -fr ${outpath}/${rf_model}
total=$(find "${fasta_path}" -name "*.toplevel.fa.gz" -not -path "*/dna_index/*" | wc -l)
find "${fasta_path}" -name "*.toplevel.fa.gz" -not -path "*/dna_index/*" -print0 | while IFS= read -r -d '' g; do
  genome_tmp=$(basename $g)
  gf=${genome_tmp%.gz}

  if [ ! -f ${outpath}/${gf}.${rf_model}.tsv ]; then
    echo "##### $COUNTER/${total}: running cmsearch on ${gf}"
    echo "unzipping ${g}"
    zcat $g > ${tmp}/${gf}

    echo "running cmsearch"
    ${cm} \
    --cpu 24 --notrunc -E 0.1 \
    -o ${outpath}/${gf}.${rf_model}.out \
    --tblout ${outpath}/${gf}.${rf_model}.tbl --fmt 3 \
    ${rf_model_path} ${tmp}/${gf}

    perl -e '
    open(my $FH, "<'${outpath}'/'${gf}'.'${rf_model}'.tbl ");
    my @lines = <$FH>;
    chomp(@lines);
    foreach my $l (@lines){
      if($l !~ /^\#/){
        my @fields = split(/\s+/,$l);
        my $str = join("\t",@fields);
        print "'$gf'\t$str\n";
      }
    }
    ' > ${outpath}/${gf}.${rf_model}.tsv
    echo "cleaning temp"
    rm ${tmp}/${gf}
  else
    echo "+++++ $COUNTER/361: results already exists: ${outpath}/${gf}.${rf_model}.tsv"
  fi

  COUNTER=$[$COUNTER +1]
done

echo "combining cmsearch"
cat ${outpath}/*.${rf_model}.tsv > ${outpath}/all.${rf_model}.tsv

echo "adding specie column"
perl -ne '
chomp($_);
my($fp,@toks) = split("\t",$_);
my($n,@rest) = split("\\.",$fp);
print "$n\t" . $_ . "\n";
' ${outpath}/all.${rf_model}.tsv > ${outpath}/all.${rf_model}.sp.tsv

echo "populate db"
sqlite3 ${outpath}/$sqlitedb '
drop table if exists rs_ensembl_genomes_'${db_rf_model}';
CREATE TABLE rs_ensembl_genomes_'${db_rf_model}'(
  specie TEXT,
  fa_path TEXT,
  target_name TEXT,
  accession TEXT,
  query_name TEXT,
  query_accession TEXT,
  model_type TEXT,
  model_from_coord TEXT,
  model_to_coord TEXT,
  target_from_coord TEXT,
  target_to_coord TEXT,
  strand TEXT,
  trunc TEXT,
  pass TEXT,
  gc TEXT,
  bias TEXT,
  score TEXT,
  evalue TEXT,
  inc TEXT,
  mdl_len TEXT,
  seq_len TEXT,
  description TEXT,
  description2 TEXT,
  description3 TEXT
);
'

sqlite3 ${outpath}/$sqlitedb '.separator "\t"' ".import ${outpath}/all.${rf_model}.sp.tsv rs_ensembl_genomes_${db_rf_model}"

sqlite3 ${outpath}/$sqlitedb "
drop index if exists rs_ensembl_genomes_${db_rf_model}_target_name_idx;
create index rs_ensembl_genomes_${db_rf_model}_idx on rs_ensembl_genomes_${db_rf_model}(description2);
"

sqlite3 ${outpath}/$sqlitedb "
alter table rs_ensembl_genomes_${db_rf_model} add column strand_lbl integer default -1;
update rs_ensembl_genomes_${db_rf_model}
set
  strand_lbl = 1
where
  strand='+';
"

sqlite3 ${outpath}/$sqlitedb "
update rs_ensembl_genomes_${db_rf_model}
set
  strand_lbl = -1
where
  strand='-';
"

echo "prepare params for overlapping feature"
python3 ${repo}/find_riboswitch_genbank_file.py \
--num_processes 4 --rf_model ${rf_model} \
--gb_path ${gb_path} \
--input_file ${outpath}/all.${rf_model}.sp.tsv \
--output_file ${outpath}/all.${rf_model}.sp.out.tsv

echo "sort by unique genomes"
sort ${outpath}/all.${rf_model}.sp.out.tsv | uniq > ${outpath}/all.${rf_model}.sp.out.tsv.uniq
awk '{print $NF,$0}' ${outpath}/all.${rf_model}.sp.out.tsv.uniq | sort | cut -f2- -d' ' \
> ${outpath}/all.${rf_model}.sp.out.tsv.uniq.genomesorted

echo "find overlapping features"
python3 ${repo}/find_riboswitch_overlapping_features_by_genome.py \
-i ${outpath}/all.${rf_model}.sp.out.tsv.uniq.genomesorted \
-o ${outpath}/all.${rf_model}.sp.out.tsv.uniq.genomesorted.over.tsv

echo "add specie column"
perl -ne '
chomp($_);
my($fp,@toks) = split("\t",$_);
my @n_tmp = split("\/",$fp);
my($n,@rest) = split("\\.",$n_tmp[-1]);
print "$n\t" . $_ . "\n";
' ${outpath}/all.${rf_model}.sp.out.tsv.uniq.genomesorted.over.tsv > ${outpath}/all.${rf_model}.sp.out.tsv.uniq.genomesorted.over.sp.tsv

echo "add ovelapping features to db"
sqlite3 ${outpath}/$sqlitedb '
drop table if exists rs_ensembl_genomes_annot_strandspec_'${db_rf_model}';
CREATE TABLE rs_ensembl_genomes_annot_strandspec_'${db_rf_model}'(
  specie TEXT,
  genome_file TEXT,
  chr TEXT,
  start INTEGER,
  end INTEGER,
  strand INTEGER,
  gene TEXT,
  loc TEXT,
  loc_strand INTEGER,
  xref TEXT,
  prod TEXT
);
'
sqlite3 ${outpath}/$sqlitedb '.separator "\t"' ".import ${outpath}/all.${rf_model}.sp.out.tsv.uniq.genomesorted.over.sp.tsv rs_ensembl_genomes_annot_strandspec_${db_rf_model}"

sqlite3 ${outpath}/$sqlitedb '
ALTER TABLE rs_ensembl_genomes_annot_strandspec_'${db_rf_model}' ADD COLUMN overlap_type TEXT DEFAULT "";
'
sqlite3 ${outpath}/$sqlitedb "
UPDATE rs_ensembl_genomes_annot_strandspec_${db_rf_model}
SET
  overlap_type=substr(loc,1,instr(loc,':')-1);
"

sqlite3 ${outpath}/$sqlitedb '.separator "\t"' '.header off' "
select distinct
  genome_file,
  chr,
  start,
  end,
  strand,
  gene
from
  rs_ensembl_genomes_annot_strandspec_${db_rf_model}
" > ${outpath}/all.${rf_model}.ensembl.genome.riboswitch.filtered.tsv

echo "get riboswitch sequences"
python ${repo}/find_riboswitch_sequence.py \
-i ${outpath}/all.${rf_model}.ensembl.genome.riboswitch.filtered.tsv \
-o ${outpath}/all.${rf_model}.ensembl.genome.riboswitch.filtered.outseq.tsv \
-g ${fasta_path}

echo "readd specie column"
perl -ne '
chomp($_);
my($fp,@toks) = split("\t",$_);
my @n_tmp = split("\/",$fp);
my($n,@rest) = split("\\.",$n_tmp[-1]);
print "$n\t" . $_ . "\n";
' ${outpath}/all.${rf_model}.ensembl.genome.riboswitch.filtered.outseq.tsv > ${outpath}/all.${rf_model}.ensembl.genome.riboswitch.filtered.outseq.sp.tsv

echo "populate db with sequence info"
sqlite3 ${outpath}/$sqlitedb '
drop table if exists rs_ensembl_genomes_annot_strandspec_seq_'${db_rf_model}';
CREATE TABLE rs_ensembl_genomes_annot_strandspec_seq_'${db_rf_model}'(
  specie TEXT,
  genome_file TEXT,
  chr TEXT,
  start INTEGER,
  end INTEGER,
  strand INTEGER,
  seq TEXT
);
'
sqlite3 ${outpath}/$sqlitedb '.separator "\t"' ".import ${outpath}/all.${rf_model}.ensembl.genome.riboswitch.filtered.outseq.sp.tsv rs_ensembl_genomes_annot_strandspec_seq_${db_rf_model}"

echo "generate report ${rf_model}.ensembl.genome.riboswitch.report.tsv"
sqlite3 ${outpath}/$sqlitedb '.separator "\t"' '.header on' "
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
  substr(ann.gene, 3, length(ann.gene) - 4) ann_gene,
  db2.g.id,
  db2.g.info,
  db2.g.gb_info,
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
  left join db2.gene_info g on db2.g.id=substr(ann.gene, 3, length(ann.gene) - 4)
  left join rs_ensembl_genomes_annot_strandspec_seq_${db_rf_model} seq
    on seq.specie=rs.specie
      AND seq.chr=rs.target_name
      AND seq.start=rs.target_from_coord
      AND seq.end=rs.target_to_coord
      AND seq.strand=rs.strand_lbl
" > ${outpath}/${rf_model}.ensembl.genome.riboswitch.report.tsv

sqlite3 ${outpath}/$sqlitedb '.separator "\t"' '.header on' "
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
  substr(ann.gene, 3, length(ann.gene) - 4) ann_gene,
  db2.g.id,
  db2.g.info,
  db2.g.gb_info,
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
  join db2.gene_info g on db2.g.id=substr(ann.gene, 3, length(ann.gene) - 4)
  left join rs_ensembl_genomes_annot_strandspec_seq_${db_rf_model} seq
    on seq.specie=rs.specie
      AND seq.chr=rs.target_name
      AND seq.start=rs.target_from_coord
      AND seq.end=rs.target_to_coord
      AND seq.strand=rs.strand_lbl
" > ${outpath}/${rf_model}.ensembl.genome.riboswitch.report.geneoverlap.tsv
