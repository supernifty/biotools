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

def get_fh(fn):
  if fn.endswith('.gz'):
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'r')

def main(bed, name, intersect, chrom_name, pos_name, padding):
  logging.info('reading %s...', bed)

  annotations = {}
  for count, line in enumerate(get_fh(bed)):
    fields = line.strip('\r\n').split('\t')
    chr = fields[0].replace('chr', '')
    if chr not in annotations:
      annotations[chr] = intervaltree.IntervalTree()
      logging.info('adding %s to annotations', chr)
    
    annotations[chr][int(fields[1])-padding:int(fields[2])+padding] = fields[3]
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
    if first: # tsv header
      first = False
      if name is None:
        sys.stdout.write('{}\n'.format(line.strip('\r\n')))
      else:
        sys.stdout.write('{}\t{}\n'.format(line.strip('\r\n'), name))
      pos = line.strip('\r\n').split('\t').index(pos_name)
      chrom = line.strip('\r\n').split('\t').index(chrom_name)
      continue

    c = line.strip('\r\n').split('\t')[chrom].replace('chr', '')
    p = line.strip('\r\n').split('\t')[pos]
    if p == pos_name: # another heading
      continue
    try:
      p = int(p)
    except:
      logging.warn('position "%s" on line %i is not numeric', p, count + 1)
      continue
    if c in annotations and len(annotations[c][p]) > 0:
      if name is None:
        sys.stdout.write('{}\n'.format(line.strip('\r\n')))
      else:
        values = set([overlap.data for overlap in annotations[c][p]])
        value = ','.join(values)
        sys.stdout.write('{}\t{}\n'.format(line.strip('\r\n'), value))
      annotated += 1
    else:
      if not intersect:
        if name:
          sys.stdout.write('{}\t{}\n'.format(line.strip('\r\n'), '-'))
        else: # this option is pointless - no name or intersect
          sys.stdout.write('{}\n'.format(line.strip('\r\n')))

    if (count + 1) % 100000 == 0:
      logging.info('%i lines %i annotations...', count + 1, annotated)

  logging.info('done. annotated %i of %i lines', annotated, count)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Add column to tsv from bed file')
  parser.add_argument('--bed', required=True, help='details to annotate')
  parser.add_argument('--padding', required=False, type=int, default=0, help='padding on bed')
  parser.add_argument('--name', required=False, help='info field name')
  parser.add_argument('--intersect', action='store_true', help='intersect as well')
  parser.add_argument('--chrom_name', required=False, default='CHROM', help='chrom column')
  parser.add_argument('--pos_name', required=False, default='POS', help='pos column')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bed, args.name, args.intersect, args.chrom_name, args.pos_name, args.padding)

