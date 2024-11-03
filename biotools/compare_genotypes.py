#!/usr/bin/env python
'''
  predict lynch pathway
'''

import argparse
import collections
import csv
import logging
import sys

def main(g1, g2, ofh):
  logging.info('starting %s...', g1)
  g1rows = []
  g1s = {}
  for r in csv.DictReader(open(g1, 'rt'), delimiter='\t'):
    g1rows.append(r)
    g1s[r['rsid']] = len(g1rows) - 1

  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=['rsid', 'g1', 'g2', 'match'])
  odw.writeheader()

  logging.info('starting %s...', g2)
  skipped = overlap = match = nomatch = 0
  for r in csv.DictReader(open(g2, 'rt'), delimiter='\t'):
    if r['rsid'] not in g1s:
      skipped += 0
    else:
      overlap += 1
      g1row = g1rows[g1s[r['rsid']]]
      is_match = r['genotype'] == g1row['genotype'] or r['genotype'] == ''.join(g1row['genotype'][::-1])
      if is_match:
        match += 1
      else:
        nomatch += 1
      odw.writerow({'rsid': r['rsid'], 'g1': g1row['genotype'], 'g2': r['genotype'], 'match': 'yes' if is_match else 'no'})
    
  logging.info('done. skipped %i overlap %i match %i nomatch %i pctmatch %.2f', skipped, overlap, match, nomatch, match / overlap)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--g1', help='more logging')
  parser.add_argument('--g2', help='more logging')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.g1, args.g2, sys.stdout)
