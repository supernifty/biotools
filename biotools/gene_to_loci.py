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

def main(genes, refseq, ignore_alt_chrom, sort, unformatted):
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
      if not ignore_alt_chrom or '_' not in row['chrom']:
        mins[row['name2']] = min([mins[row['name2']], int(row['txStart'])])
        maxs[row['name2']] = max([maxs[row['name2']], int(row['txEnd'])])
        chroms[row['name2']] = row['chrom']
      else:
        logging.debug('skipped alternative chrom %s', row['chrom'])

  if unformatted:
    # assume sort
    items = []
    for gene in genes:
      if gene in chroms:
        items.append((gene, chroms[gene], mins[gene], maxs[gene]))
    sys.stdout.write(' '.join(['{}:{}-{}'.format(item[1], item[2], item[3]) for item in sorted(items, key=lambda item: (int(item[1].replace('chr', '')), item[2]))]))
  else:
    sys.stdout.write('Gene\tLoci\n')
    if sort:
      items = []
      for gene in genes:
        if gene in chroms:
          items.append((gene, chroms[gene], mins[gene], maxs[gene]))
      for item in sorted(items, key=lambda item: (int(item[1].replace('chr', '')), item[2])):
        sys.stdout.write('{}\t{}:{}-{}\n'.format(item[0], item[1], item[2], item[3]))
  
    else:
      for gene in sorted(genes):
        if gene not in chroms:
          logging.warn('Failed to find %s in refseq', gene)
        else:
          sys.stdout.write('{}\t{}:{}-{}\n'.format(gene, chroms[gene], mins[gene], maxs[gene]))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate region for specified genes')
  parser.add_argument('--genes', required=True, nargs='+', help='')
  parser.add_argument('--refseq', required=True, help='refseq gz')
  parser.add_argument('--ignore_alt_chroms', action='store_true', help='refseq gz')
  parser.add_argument('--sort', action='store_true', help='sort on chromosome then position')
  parser.add_argument('--unformatted', action='store_true', help='write as space separated list of regions')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genes, args.refseq, args.ignore_alt_chroms, args.sort, args.unformatted)
