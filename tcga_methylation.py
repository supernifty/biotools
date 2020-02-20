#!/usr/bin/env python

import argparse
import collections
import logging
import sys

MLH1_PROBES=set(['cg23658326', 'cg11600697', 'cg21490561'])

#CACNA1G = c("cg23614229", "cg11262815")
#RUNX3 = c("cg19590532", "cg06377275"),
#SOCS1 = "cg06220235"
#NEUROG1 = c("cg07781035","cg04620091", "cg22630755")
#IGF2 = "cg16977706"

CIMP_CUTOFF=0.2
CIMP_PROBES= {
  'CACNA1G': set(["cg23614229", "cg11262815"]),
  'RUNX3': set(["cg19590532", "cg06377275"]),
  'SOCS1': set(["cg06220235"]),
  'NEUROG1': set(["cg07781035","cg04620091", "cg22630755"]),
  'IGF2': set(["cg16977706"])
}

def main(files):
  sys.stdout.write('Sample\tMLH1\tCIMP\t{}\t{}\n'.format('\t'.join(sorted(MLH1_PROBES)), '\t'.join(['\t'.join(['{}_{}'.format(gene, x) for x in sorted(CIMP_PROBES[gene])]) for gene in sorted(CIMP_PROBES)])))
  for f in open(files, 'r'):
    f = f.strip('\n')
    logging.info('processing %s...', f)
    mlh1_values = {}
    cimp_values = collections.defaultdict(list)
    probes = {}
    for line in open(f, 'r'):
      fields = line.strip('\n').split('\t')
      if fields[0] in MLH1_PROBES:
        mlh1_values[fields[0]] = fields[1]
      for gene in CIMP_PROBES:
        if fields[0] in CIMP_PROBES[gene]:
          cimp_values[gene].append(fields[1])
          probes[fields[0]] = fields[1]

    sample = 'TCGA{}'.format(f.split('TCGA')[-1].replace('.gdc_hg38.txt', ''))
    sys.stdout.write('{}\t{:.2f}\t'.format(sample, sum([float(mlh1_values[x]) for x in mlh1_values]) / len(mlh1_values)))
    cimp_count = 0
    for gene in CIMP_PROBES:
      if sum([float(x) for x in cimp_values[gene]]) / len(cimp_values[gene]) > CIMP_CUTOFF:
        cimp_count += 1

    sys.stdout.write('{}\t'.format(cimp_count))

    # mlh1 probes
    sys.stdout.write('\t'.join(mlh1_values.get(x, 'NA') for x in sorted(MLH1_PROBES)))

    # cimp probes
    for gene in sorted(CIMP_PROBES):
      for probe in sorted(CIMP_PROBES[gene]):
        sys.stdout.write('\t{}'.format(probes.get(probe, 'NA')))
 
    sys.stdout.write('\n')

    logging.info('processing %s: done', f)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--files', required=True, help='file containing filenames')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.files)

