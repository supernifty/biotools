#!/usr/bin/env python
'''
  take fields from vcf
'''

import argparse
import collections
import logging
import os.path
import sys

import cyvcf2

C=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']


def main(vcfs, ofh, indels_only):

  v = {}
  for kv in vcfs:
    name, vcf = kv.split('=')
    logging.info('reading %s with name %s...', vcf, name)

    count = 0
    for count, variant in enumerate(cyvcf2.VCF(vcf)):
      if indels_only and len(variant.REF) == 1 and len(variant.ALT) == 1 and len(variant.ALT[0]) == 1:
        continue
      if variant.CHROM not in v:
        v[variant.CHROM] = collections.defaultdict(list)
        logging.info('adding %s to annotations', variant.CHROM)
      v[variant.CHROM][(variant.POS, variant.REF, ','.join(variant.ALT))].append({'caller': name})
      if (count + 1) % 100000 == 0:
        logging.info('%i lines...', count + 1)

    logging.info('reading %s: done with %i', vcf, count)

  # write result
  ofh.write('''##fileformat=VCFv4.1
##INFO=<ID=caller,Number=1,Type=String,Description="caller(s)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
''')

  written = 0
  for c in C:
    if c not in v:
      continue
    logging.info('writing %s...', c)
    for k in sorted(v[c].keys()):
      logging.debug('writing key %s', k)
      ofh.write('{}\t{}\t.\t{}\t{}\t.\tPASS\tcaller={}\n'.format(c, k[0], k[1], k[2], ','.join([x['caller'] for x in v[c][k]])))
      written += 1

  logging.info('done writing %i', written)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Merge two vcfs')
  parser.add_argument('--vcfs', required=True, nargs='+', help='vcf files to merge of the form name=filename')
  parser.add_argument('--indels_only', action='store_true', help='indels only')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, sys.stdout, args.indels_only)

