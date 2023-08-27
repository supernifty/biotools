#!/usr/bin/env python
'''
  clinvar stats
'''

import argparse
import collections
import logging
import sys

import cyvcf2

def main(gene_filter):
  logging.info('starting...')

  # 1       45794974        479982  G       C       .       .       AF_EXAC=0.00001;ALLELEID=472328;CLNDISDB=MedGen:C0027672,SNOMED_CT:699346009|MedGen:CN517202;CLNDN=Hereditary_cancer-predisposing_syndrome|not_provided;CLNHGVS=NC_000001.10:g.45794974G>C;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=MUTYH:4595;MC=SO:0001619|non-coding_transcript_variant,SO:0001624|3_prime_UTR_variant;ORIGIN=1;RS=758118037
  stat = collections.defaultdict(int)
  genes = collections.defaultdict(int)

  count = skipped = 0
  for count, variant in enumerate(cyvcf2.VCF('-')):
    try:
      gene = variant.INFO['GENEINFO']
    except KeyError:
      gene = 'notspecified'

    genes[gene] += 1

    gene_name = gene.split(':')[0]

    if gene_filter is not None and gene_name != gene_filter:
      continue 
    try:
      stat[variant.INFO['CLNSIG']] += 1
    except:
      # no clnsig
      stat['missing'] += 1
      skipped += 1
    if count % 10000 == 0:
      logging.debug('%s processed...', count)

  sys.stdout.write('CLNSIG\tCount\tPct\n')
  total = sum([stat[c] for c in stat])
  for c in stat:
    sys.stdout.write('{}\t{}\t{:.3f}\n'.format(c, stat[c], stat[c] / total))

  #logging.info(genes)
  logging.info('done processing %i skipped %i', count, skipped)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='clinvar stats')
  parser.add_argument('--gene', required=False, help='stats for this gene')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.gene)

