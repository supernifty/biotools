#!/usr/bin/env python
'''
  take fields from vcf
'''

import argparse
import csv
import collections
import gzip
import logging
import sys

import cyvcf2
import intervaltree

def main(refseq):
  logging.info('reading %s...', refseq)
  tree = {}
  # #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
  #0       NM_001276352.2  chr1    -       67092164        67134970        67093579        67127240        9       67092164,67096251,67103237,67111576,67115351,67125751,67127165,67131141,67134929,       67093604,67096321,67103382,67111644,67115464,67125909,67127257,67131227,67134970, 0       C1orf141        cmpl    cmpl    2,1,0,1,2,0,0,-1,-1,
  for row in csv.DictReader(gzip.open(refseq, 'rt'), delimiter='\t'):
    c = row['chrom']
    g = row['name2']
    s = row['txStart']
    e = row['txEnd']
    if c not in tree:
      tree[c] = intervaltree.IntervalTree()
    tree[c][int(s):int(e)] = g
  logging.info('reading %s: done', refseq)

  # now annotate input
  logging.info('reading from stdin...')
  vcf_in = cyvcf2.VCF('-')
  vcf_in.add_info_to_header({'ID': 'gene', 'Description': 'gene according to refseq', 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)

  count = 0
  for count, variant in enumerate(vcf_in):
    if variant.CHROM in tree:
      genes = ','.join(sorted(list(set([overlap.data for overlap in tree[variant.CHROM].at(variant.POS)]))))
      variant.INFO['gene'] = genes
    sys.stdout.write(str(variant))
    if (count + 1) % 10000 == 0:
      logging.info('reading %s: %i lines processed', variant.CHROM, count + 1)
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Add gene from refseq to vcf')
  parser.add_argument('--refseq', required=True, help='refseq file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.refseq)

