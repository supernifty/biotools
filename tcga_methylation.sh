#!/bin/bash
#SBATCH --job-name=team
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=16:00:00
#SBATCH --partition=physical,snowy
#SBATCH --account=punim0567
#SBATCH --mem=2G

set -o errexit

OUT="methylation.tsv"
ROOT=/data/cephfs/punim0567/data/public_datasets/tcga/
#/data/cephfs/punim0567/data/public_datasets/tcga/./be76b382-ebff-4349-b55a-4b85cd98e468/jhu-usc.edu_COAD.HumanMethylation450.2.lvl-3.TCGA-CM-6161-01A-11D-1651-05.gdc_hg38.txt
echo "finding files..."
FILES=methylation.coadread.txt
grep "jhu-usc.edu_COAD.HumanMethylation450.*.txt$" $ROOT/files | sed 's/^/\/data\/cephfs\/punim0567\/data\/public_datasets\/tcga\//' > $FILES
grep "jhu-usc.edu_READ.HumanMethylation450.*.txt$" $ROOT/files | sed 's/^/\/data\/cephfs\/punim0567\/data\/public_datasets\/tcga\//' >> $FILES
echo "finding files: done"

python /data/cephfs/punim0567/peter/src/biotools/tcga_methylation.py --files $FILES > $OUT

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
