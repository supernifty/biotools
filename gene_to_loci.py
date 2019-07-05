#!/usr/bin/env python
'''
  list of genes to loci given refseq
  e.g.
  python gene_to_loci.py --genes MSH2 RNF43 --refseq ~/crc/data/public_datasets/refseq.ucsc.hg19.180829.gz
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

def main(genes, refseq):
  logging.info('starting...')
  genes = set(genes)
  mins = collections.defaultdict(int)
  mins.default_factory = lambda: 3e9
  maxs = collections.defaultdict(int)
  chroms = {}
  for row in csv.DictReader(gzip.open(refseq, 'rt'), delimiter='\t'):
    if row['name2'] in genes:
      if row['name2'] in chroms and row['chrom'] != chroms[row['name2']]:
        logging.warn('multiple chromosomes found for %s', row['name2'])
      mins[row['name2']] = min([mins[row['name2']], int(row['txStart'])])
      maxs[row['name2']] = max([maxs[row['name2']], int(row['txEnd'])])
      chroms[row['name2']] = row['chrom']

  sys.stdout.write('Gene\tLoci\n')
  for gene in genes:
    if gene not in chroms:
      logging.warn('Failed to find %s in refseq', gene)
    else:
      sys.stdout.write('{}\t{}:{}-{}\n'.format(gene, chroms[gene], mins[gene], maxs[gene]))
    
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate region for specified genes')
  parser.add_argument('--genes', required=True, nargs='+', help='')
  parser.add_argument('--refseq', required=True, help='refseq gz')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genes, args.refseq)
