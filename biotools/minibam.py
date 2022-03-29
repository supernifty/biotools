#!/usr/bin/env python
'''
  generate a minibam based on a vcf
  requires samtools
  e.g.
  python minibam.py --bam bamfile --out minibam < vcf
'''

import argparse
import collections
import csv
import gzip
import logging
import os
import sys
import tempfile

import cyvcf2

def main(bam, minibam, padding):
  logging.info('reading from stdin...')
  items = []

  fp = tempfile.NamedTemporaryFile(mode='w', delete=False)
  c = 0
  for c, r in enumerate(cyvcf2.VCF('-')):
    if '_' in r.CHROM: # skip alt chroms
      continue
    #items.append((r.CHROM, r.POS-padding, r.POS+padding))
    fp.write('{}\t{}\t{}\n'.format(r.CHROM, r.POS-padding, r.POS+padding))
  fp.close()
  
  logging.info('%s variants', c)

  # samtools view -h -b {input.bam} {config[regions_of_interest]} > {output}
  #cmd = 'samtools view -h -b {} -L {} > {}'.format(bam, fp.name, minibam)
  cmd = 'bedtools intersect -a {} -b {} > {}'.format(bam, fp.name, minibam)
  logging.debug('command is: %s', cmd)
  if os.system(cmd) == 0:
    logging.info('success!')
  else:
    logging.warn('command failed')

  os.remove(fp.name)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate region from specified VCF')
  parser.add_argument('--padding', required=False, default=100, help='padding around variants')
  parser.add_argument('--bam', required=True, help='source bam')
  parser.add_argument('--out', required=True, help='minbam out')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bam, args.out, args.padding)
