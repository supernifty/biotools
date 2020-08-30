#!/bin/bash
#SBATCH --job-name=methylation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --partition=physical,snowy
#SBATCH --account=punim0567
#SBATCH --mem=2G

set -o errexit

# have been running this in /data/gpfs/projects/punim0567/peter/analyses/mutational_signature_search_review
cd /data/gpfs/projects/punim0567/peter/analyses/mutational_signature_search_review

OUT_450="methylation.450.tsv"
OUT_27="methylation.27.tsv"
OUT="methylation.tsv"
ROOT=/data/cephfs/punim0567/data/public_datasets/tcga/

# 450k
#/data/cephfs/punim0567/data/public_datasets/tcga/./be76b382-ebff-4349-b55a-4b85cd98e468/jhu-usc.edu_COAD.HumanMethylation450.2.lvl-3.TCGA-CM-6161-01A-11D-1651-05.gdc_hg38.txt

# 27 
# /data/cephfs/punim0567/data/public_datasets/tcga/./3c2dae10-87fe-46f6-941a-b8be6e024b5b/jhu-usc.edu_COAD.HumanMethylation27.7.lvl-3.TCGA-AY-4070-01A-01D-1110-05.gdc_hg38.txt
# /data/cephfs/punim0567/data/public_datasets/tcga/./019a5ba8-0368-4904-a628-a9373d838c09/jhu-usc.edu_COAD.HumanMethylation27.5.lvl-3.TCGA-AA-3680-01A-01D-0904-05.gdc_hg38.txt

echo "$(date) finding files..."
FILES_450k=methylation.coadread.450k.txt
FILES_27k=methylation.coadread.27k.txt
grep "jhu-usc.edu_COAD.HumanMethylation450.*.txt$" $ROOT/files | sed 's/^/\/data\/cephfs\/punim0567\/data\/public_datasets\/tcga\//' > $FILES_450k
grep "jhu-usc.edu_READ.HumanMethylation450.*.txt$" $ROOT/files | sed 's/^/\/data\/cephfs\/punim0567\/data\/public_datasets\/tcga\//' >> $FILES_450k
grep "jhu-usc.edu_COAD.HumanMethylation27.*.txt$" $ROOT/files | sed 's/^/\/data\/cephfs\/punim0567\/data\/public_datasets\/tcga\//' > $FILES_27k
grep "jhu-usc.edu_READ.HumanMethylation27.*.txt$" $ROOT/files | sed 's/^/\/data\/cephfs\/punim0567\/data\/public_datasets\/tcga\//' >> $FILES_27k
echo "$(date) finding files: done"

#python /data/cephfs/punim0567/peter/src/biotools/tcga_methylation.py --files $FILES_27k --array 27k --verbose > $OUT_27
#python /data/cephfs/punim0567/peter/src/biotools/tcga_methylation.py --files $FILES_450k --array 450k > $OUT_450

csvcols.py --delimiter '	' --cols Sample MLH1 CIMP < $OUT_450 | csvadd.py --delimiter '	' --name Array --value 450k | \
  sed '/^....-..-....-[^0]/d' > tmp.methylation # exclude non-tumours
csvcols.py --delimiter '	' --cols Sample MLH1 CIMP < $OUT_27 | csvadd.py --delimiter '	' --name Array --value 27k | sed '1d' | \
  sed '/^....-..-....-[^0]/d' >> tmp.methylation # exclude non-tumours

csvsort.py --delimiter '	' --cols Sample < tmp.methylation |\
  csvadd.py --delimiter '	' --name M1 --value No --rule 'MLH1>0.2:Yes' |\
  csvadd.py --delimiter '	' --name M2 --value No --rule 'CIMP>2:Yes' |\
  csvop.py --delimiter '	' --dest M1M2 --op concat --cols M1 M2 --join_string ',' |\
  csvadd.py --delimiter '	' --name Outcome --value No --rule M1M2=Yes,Yes:Yes > $OUT

# want to only include 27k that is not already provided as 450k
python src/check_meth.py < $OUT 

echo "$(date) calculation: done"

#>$OUT
#echo "Sample	
#for g in $(cat $FILES); do
#  f=$ROOT$g
#  echo "$f..."
#  result=$(python tcga_methylation.py < $f)
#  echo "$f: done"
#  s=${f/*TCGA/TCGA}
#  t=${s/.gdc_hg38.txt/}
#  echo "$t	$result" >> $OUT
#done

cd -
