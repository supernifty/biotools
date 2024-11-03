# biotools
Python based utility scripts

## annotate vcf

```
python annotate_vcf.py --vcf clinvar_20191219.vcf.gz --fields CLNDN CLNSIG < in.vcf.gz | bgzip > out.vcf.gz
```

## compare_genotypes.py
find percentage overlap between two genotype files

```
python src/compare_genotypes.py --g1 out/23andMe.genotype --g2 out/p.vcf.genotype > out/p.comparison
```

## convert_vcf_to_genotype.py
generates a genotype file given a vcf

```
python biotools/convert_vcf_to_genotype.py --vcf in/p.vcf.gz > out/p.vcf.genotype
```

## find variants

## count indexes in fastq files

