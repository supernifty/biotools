#!/usr/bin/env python
'''
  list of genes to loci given refseq
  e.g.
  python vcf_to_loci.py --padding 100 < vcf
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

import cyvcf2

def main(padding):
  logging.info('starting...')
  items = []

  for r in cyvcf2.VCF('-'):
    if '_' in r.CHROM: # skip alt chroms
      continue
    items.append((r.CHROM, r.POS-padding, r.POS+padding))

  result = []
  last = ('', 0)
  for item in sorted(items, key=lambda item: (int(item[0].replace('chr', '')), item[1])):
    if item[0] == last[0] and item[1] < last[1]:
      continue
    result.append('{}:{}-{}'.format(item[0], item[1], item[2]))
    last = (item[0], item[2])

  sys.stdout.write(' '.join(result))
  sys.stdout.write('\n')
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate region from specified VCF')
  parser.add_argument('--padding', required=False, default=100, help='padding around variants')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.padding)
