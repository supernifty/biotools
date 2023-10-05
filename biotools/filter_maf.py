#!/usr/bin/env python
'''
  filter on bed
'''

import argparse
import csv
import logging
import sys

import intervaltree
import cyvcf2

def bed_to_tree(bed):
  logging.info('parsing {}...'.format(bed))
  tree = {}
  size = overlaps = skipped = included = 0
  if bed.endswith('.gz'):
    fh = gzip.open(bed, 'rt')
  else:
    fh = open(bed, 'r')
  for line_count, line in enumerate(fh):
    fields = line.strip('\n').split('\t')
    if len(fields) < 3:
      skipped += 1
      continue
    chrom, start, finish = fields[:3]
    if chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
      logging.info('processing chromosome %s...', chrom)
    s = int(start)
    f = int(finish)
    # does it overlap with itself?
    overlap = tree[chrom].overlap(s, f)
    if len(overlap) == 0:
      size += f - s
      tree[chrom][s:f] = True
      included += 1
    else:
      if len(overlap) > 1:
        logging.warn('unexpected multiple overlaps: %s', overlap)
      for item in overlap:
        #logging.debug('overlap %s:%i-%i with %s:%i-%i', chrom, item.begin, item.end, chrom, s, f)
        new_begin = min(s, item.begin)
        new_end = max(f, item.end)
        tree[chrom].remove(item)
        overlaps += 1

      tree[chrom][new_begin:new_end] = True
      size += (new_end - new_begin) - (f - s)
      s = new_begin
      f = new_end

    if line_count % 100000 == 0:
      logging.debug('parsing {}: {} lines parsed. skipped {}. {} overlaps. size {}'.format(bed, line_count, skipped, overlaps, size))
  logging.info('parsing {}: done. lines skipped: {}. size: {}. count: {}'.format(bed, skipped, size, included))
  return tree, size

def main(ifh, ofh, bed):
  logging.info('reading %s...', bed)
  tree, size = bed_to_tree(bed)

  logging.info('reading maf from stdin...')
  idr = csv.DictReader(ifh, delimiter='\t')
  odw = csv.DictWriter(ofh, delimiter='\t', fieldnames=idr.fieldnames)
  odw.writeheader()
  included = total = 0
  for r in idr:
    total += 1
    # Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position
    chrom = r['Chromosome']
    if chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom not in tree:
      continue
    p = int(r['Start_Position'])
    if len(tree[chrom].at(p-1)) > 0:
      odw.writerow(r)
      included += 1

  logging.info('included %i of %i', included, total)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter MAF on bed')
  parser.add_argument('--bed', required=True, help='bed to filter')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  # sample af vcf
  main(sys.stdin, sys.stdout, args.bed)
