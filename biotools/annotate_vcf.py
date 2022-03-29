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

def main(vcf, fields):
  logging.info('reading %s...', vcf)

  annotations = {}
  for count, variant in enumerate(cyvcf2.VCF(vcf)):
    chr = variant.CHROM.replace('chr', '')
    if chr not in annotations:
      annotations[chr] = {}
      logging.info('adding %s to annotations', chr)
    try:
      annotations[chr]['{}/{}/{}'.format(variant.POS, variant.REF, variant.ALT)] = {name: variant.INFO[name] for name in fields}
    except KeyError:
      logging.debug('line %i %s:%i: fields not found', count, chr, variant.POS)
    if (count + 1) % 100000 == 0:
      logging.info('%i lines...', count + 1)
  logging.debug('reading %s: %i lines processed', chr, count + 1)

  logging.info('reading %s: done', vcf)

  # now annotate input
  logging.info('reading from stdin...')

  vcf_in = cyvcf2.VCF('-')
  for field in fields:
    vcf_in.add_info_to_header({'ID': field, 'Description': 'Annotated field {}'.format(field), 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)

  annotated = 0
  count = 0
  seen = set()
  for count, variant in enumerate(vcf_in):
    chr = variant.CHROM.replace('chr', '')
    if chr in annotations:
      key = '{}/{}/{}'.format(variant.POS, variant.REF, variant.ALT)
      if key in annotations[chr]:
        for field in fields:
          variant.INFO[field] = annotations[chr][key][field]
        annotated += 1
    else:
      if chr not in seen:
        seen.add(chr)
        logging.warn('chromosome %s not seen in annotations', chr)
    sys.stdout.write(str(variant))
    if (count + 1) % 10000 == 0:
      logging.info('reading %s: %i lines processed %i annotated', chr, count + 1, annotated)
  
  logging.info('done. annotated %i of %i variants', annotated, count)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcf', required=True, help='details to annotate')
  parser.add_argument('--fields', required=True, nargs='+', help='info fields')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcf, args.fields)

