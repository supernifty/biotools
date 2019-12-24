#!/usr/bin/env python
'''
  take fields from vcf
'''

import argparse
import collections
import gzip
import logging
import os.path
import sys

import intervaltree

def main(bed, name):
  logging.info('reading %s...', bed)

  annotations = {}
  for count, line in enumerate(gzip.open(bed, 'rt')):
    fields = line.strip('\r\n').split('\t')
    chr = fields[0].replace('chr', '')
    if chr not in annotations:
      annotations[chr] = intervaltree.IntervalTree()
      logging.info('adding %s to annotations', chr)
    
    annotations[chr][int(fields[1]):int(fields[2])] = fields[3]
    if (count + 1) % 10000 == 0:
      logging.info('%i lines last annotation %s...', count + 1, fields[3])
  logging.debug('reading %s: %i lines processed', chr, count + 1)

  logging.info('reading %s: done', bed)

  # now annotate input
  logging.info('reading from stdin...')

  first = True
  pos = -1
  chrom = -1
  annotated = 0
  count = 0
  for count, line in enumerate(sys.stdin):
    if first:
      first = False
      sys.stdout.write('{}\t{}\n'.format(line.strip('\r\n'), name))
      pos = line.strip('\r\n').split('\t').index('POS')
      chrom = line.strip('\r\n').split('\t').index('CHROM')
      continue

    c = line.strip('\r\n').split('\t')[chrom].replace('chr', '')
    p = line.strip('\r\n').split('\t')[pos]
    if p == 'POS': # another heading
      continue
    p = int(p)
    if c in annotations and len(annotations[c][p]) > 0:
      values = set([overlap.data for overlap in annotations[c][p]])
      value = ','.join(values)
      sys.stdout.write('{}\t{}\n'.format(line.strip('\r\n'), value))
      annotated += 1
    else:
      sys.stdout.write('{}\t{}\n'.format(line.strip('\r\n'), '-'))

    if (count + 1) % 100000 == 0:
      logging.info('%i lines %i annotations...', count + 1, annotated)

  logging.info('done. annotated %i of %i lines', annotated, count)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--bed', required=True, help='details to annotate')
  parser.add_argument('--name', required=True, help='info field name')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bed, args.name)

