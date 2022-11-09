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

def skip_header(fn, zipped):
  if zipped:
    for line in gzip.open(fn, 'rt'):
      if not line.startswith('#'):
        yield line
  else:
    for line in open(fn, 'rt'):
      if not line.startswith('#'):
        yield line

def main(vcfs, mafs, variants_fn, padding, maf_chromosome, maf_pos, maf_ref, maf_alt, maf_sample, maf_zipped, count_snvs):
  logging.info('reading %s...', variants_fn)

  # variants to find
  variants = {}
  names = {}
  for r in csv.DictReader(open(variants_fn, 'rt'), delimiter='\t'): # expect chrom pos ref alt
    key = (r['chrom'], r['pos'], r['ref'], r['alt'])
    variants[key] = 0
    names[key] = r.get('name', '{}-{}-{}-{}'.format(r['chrom'], r['pos'], r['ref'], r['alt']))
    if padding > 0:
      for p in range(padding + 1):
        newpos = int(r['pos']) + p
        key = (r['chrom'], str(newpos), r['ref'], r['alt'])
        names[key] = r.get('name', '{}-{}-{}-{}'.format(r['chrom'], str(newpos), r['ref'], r['alt']))
        variants[key] = p
        newpos = int(r['pos']) - p
        key = (r['chrom'], str(newpos), r['ref'], r['alt'])
        names[key] = r.get('name', '{}-{}-{}-{}'.format(r['chrom'], str(newpos), r['ref'], r['alt']))
        variants[key] = -p

  logging.info('checking for %i variants in total', len(variants))
  logging.debug('names: %s', names)

  found = collections.defaultdict(set)

  # each vcf
  for idx, vcf in enumerate(vcfs):
    logging.info('reading %i of %i: %s...', idx + 1, len(vcfs), vcf)
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
  snvs = collections.defaultdict(int)
  all_mafs = set()
  for maf in mafs:
    logging.info('reading %s...', maf)
    for idx, row in enumerate(csv.DictReader(skip_header(maf, maf_zipped), delimiter='\t')):
      key = (row[maf_chromosome], row[maf_pos], row[maf_ref], row[maf_alt])
      where = '{}-{}'.format(maf, row[maf_sample])
      all_mafs.add(where)
      if key in variants:
        found[where].add(key)
        found_count += 1
        logging.debug('found %s in %s', key, where)
      if len(row[maf_ref]) == 1 and len(row[maf_alt]) == 1:
        snvs[where] += 1
      if idx % 10000 == 0:
        logging.debug('processed %i variants found %i last key %s', idx, found_count, key)

  logging.info('writing to stdout...')
  fieldnames = ['Filename']
  if count_snvs:
    fieldnames.append('snvs')
  ofh = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=fieldnames + [names.get(v, '/'.join(list(v))) for v in sorted(variants)])
  ofh.writeheader()
  for vcf in vcfs:
    r = {'Filename': vcf}
    for v in variants:
      k = '/'.join(list(v))
      if v in found[vcf]:
        r[names.get(v, k)] = 'Present'
      else:
        r[names.get(v, k)] = 'Absent'
    ofh.writerow(r)

  for maf in all_mafs:
    r = {'Filename': maf}
    if count_snvs:
      r['snvs'] = snvs[maf]
    for v in variants:
      k = '/'.join(list(v))
      if v in found[maf]:
        r[names.get(v, k)] = 'Present'
      else:
        r[names.get(v, k)] = 'Absent'
    ofh.writerow(r)

  # write a summary to stderr
  ofh = csv.DictWriter(sys.stderr, delimiter='\t', fieldnames=['Variant', 'Seen', 'Samples'])
  ofh.writeheader()
  samples = len(vcfs) + len(all_mafs)
  for v in sorted(variants):
    k = '/'.join(list(v))
    seen = sum([1 if v in found[x] else 0 for x in found])
    ofh.writerow({'Variant': names.get(v, k), 'Seen': seen, 'Samples': samples})
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='find variants')
  parser.add_argument('--vcfs', required=False, nargs='*', default=[], help='vcfs to check')
  parser.add_argument('--mafs', required=False, nargs='*', default=[], help='mafs to check')
  parser.add_argument('--maf_chromosome', required=False, default='Chromosome', help='maf column name for chromosome')
  parser.add_argument('--maf_pos', required=False, default='Start_Position', help='maf column name for start pos')
  parser.add_argument('--maf_ref', required=False, default='Reference_Allele', help='maf column name for ref allele')
  parser.add_argument('--maf_alt', required=False, default='Tumor_Seq_Allele2', help='maf column name for alt allele')
  parser.add_argument('--maf_sample', required=False, default='Chromosome', help='maf column name for sample')
  parser.add_argument('--maf_zipped', action='store_true', help='is the maf zipped?')
  parser.add_argument('--padding', required=False, default=1, type=int, help='file containing variants to look for')
  parser.add_argument('--variants', required=True, help='file containing variants to look for - tsv with fields chrom pos ref alt (optional name)')
  parser.add_argument('--count_snvs', action='store_true', help='more logging')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.mafs, args.variants, args.padding, args.maf_chromosome, args.maf_pos, args.maf_ref, args.maf_alt, args.maf_sample, args.maf_zipped, args.count_snvs)
