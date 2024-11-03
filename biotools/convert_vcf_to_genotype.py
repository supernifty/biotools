#!/usr/bin/env python
'''
  vcf -> 23andme
'''

import argparse
import collections
import csv
import logging
import sys

import cyvcf2

def main(vcf_fn, ofh):
  logging.info('starting...')
  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=['rsid', 'chromosome', 'position', 'genotype', 'genotype_orig'])
  odw.writeheader()
  vcf = cyvcf2.VCF(vcf_fn)
  for i, v in enumerate(vcf):
    gt = v.genotypes[0] #v.format('GT')[0]
    row = {'rsid': v.ID, 'chromosome': v.CHROM, 'position': v.POS, 'genotype_orig': gt}
    #logging.info(gt)
    row['genotype'] = '{}{}'.format(v.REF if gt[0] == 0 else v.ALT[0], v.REF if gt[1] == 0 else v.ALT[0]) 
    odw.writerow(row)
    if i % 100000 == 0:
      logging.info('%i variants', i)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcf', help='vcf')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcf, sys.stdout)
