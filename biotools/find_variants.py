#!/usr/bin/env python
'''
  make a table checking for specific variants
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

import cyvcf2

def skip_header(fn):
  for line in gzip.open(fn, 'rt'):
    if not line.startswith('#'):
      yield line

def main(vcfs, mafs, variants_fn, padding):
  logging.info('reading %s...', variants_fn)

  # variants to find
  variants = {}
  for r in csv.DictReader(open(variants_fn, 'rt'), delimiter='\t'): # expect chrom pos ref alt
    key = (r['chrom'], r['pos'], r['ref'], r['alt'])
    variants[key] = 0
    if padding > 0:
      for p in range(padding + 1):
        newpos = int(r['pos']) + p
        key = (r['chrom'], str(newpos), r['ref'], r['alt'])
        variants[key] = p
        newpos = int(r['pos']) - p
        key = (r['chrom'], str(newpos), r['ref'], r['alt'])
        variants[key] = -p

  logging.info('checking for %i variants in total', len(variants))

  found = collections.defaultdict(set)

  # each vcf
  for vcf in vcfs:
    logging.info('reading %s...', vcf)
    last_chr = None
    for variant in cyvcf2.VCF(vcf):
      key = (variant.CHROM, str(variant.POS), variant.REF, variant.ALT[0])
      if key in variants:
        found[vcf].add(key)
        logging.debug('found %s in %s', key, vcf)
      if variant.CHROM != last_chr:
        last_chr = variant.CHROM
        logging.debug('chromosome %s', last_chr)

  # Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1
  found_count = 0
  all_mafs = set()
  for maf in mafs:
    logging.info('reading %s...', maf)
    for idx, row in enumerate(csv.DictReader(skip_header(maf), delimiter='\t')):
      key = (row['Chromosome'], row['Start_Position'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
      where = '{}-{}'.format(maf, row['Tumor_Sample_Barcode'])
      all_mafs.add(where)
      if key in variants:
        found[where].add(key)
        found_count += 1
        logging.debug('found %s in %s', key, where)
      if idx % 10000 == 0:
        logging.debug('processed %i variants found %i last key %s', idx, found_count, key)

  logging.info('writing to stdout...')
  ofh = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=['Filename'] + ['/'.join(list(v)) for v in sorted(variants)])
  ofh.writeheader()
  for vcf in vcfs:
    r = {'Filename': vcf}
    for v in variants:
      k = '/'.join(list(v))
      if v in found[vcf]:
        r[k] = 'Present'
      else:
        r[k] = 'Absent'
    ofh.writerow(r)

  for maf in all_mafs:
    r = {'Filename': maf}
    for v in variants:
      k = '/'.join(list(v))
      if v in found[maf]:
        r[k] = 'Present'
      else:
        r[k] = 'Absent'
    ofh.writerow(r)
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='find variants')
  parser.add_argument('--vcfs', required=False, nargs='*', default=[], help='vcfs to check')
  parser.add_argument('--mafs', required=False, nargs='*', default=[], help='mafs to check')
  parser.add_argument('--padding', required=False, default=1, type=int, help='tolerance for position')
  parser.add_argument('--variants', required=True, help='file containing variants to look for - tsv with chrom pos ref alt')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.mafs, args.variants, args.padding)
