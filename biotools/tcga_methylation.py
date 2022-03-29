#!/usr/bin/env python

import argparse
import collections
import logging
import sys

# 450k probes
MLH1_PROBES_450k=set(['cg23658326', 'cg11600697', 'cg21490561'])

CIMP_CUTOFF_450k=0.2
CIMP_PROBES_450k= {
  'CACNA1G': set(["cg23614229", "cg11262815"]),
  'RUNX3': set(["cg19590532", "cg06377278"]),
  'SOCS1': set(["cg06220235"]),
  'NEUROG1': set(["cg07781035","cg04620091", "cg22630755"]),
  'IGF2': set(["cg16977706"])
}

#MLH1_PROBES_27k=set(['cg10990993', 'cg02279071', 'cg00893636'])
MLH1_PROBES_27k=set(['cg00893636'])

CIMP_CUTOFF_27k=0.2
CIMP_PROBES_27k= {
  'CACNA1G': set(['cg11262815', 'cg18337803', 'cg18454685']),
  'RUNX3': set(['cg00117172', 'cg06377278']),
  'SOCS1': set(['cg06220235']),
  'NEUROG1': set(['cg04897683', 'cg11248413', 'cg04330449']),
  'IGF2': set(['cg02166532', 'cg20339650'])
}

# useful script: csvcols.py --delimiter '  ' --cols 'Composite Element REF' 'Start' 'Gene_Symbol' < /data/cephfs/punim0567/data/public_datasets/tcga/./3c2dae10-87fe-46f6-941a-b8be6e024b5b/jhu-usc.edu_COAD.HumanMethylation27.7.lvl-3.TCGA-AY-4070-01A-01D-1110-05.gdc_hg38.txt | csvfilter.py --delimiter '   ' --filter 'Gene_Symbol%IGF2' | csvsort.py --delimiter '        ' --cols Start --numeric

def main(files, array):
  if array == '450k':
    mlh1_probes, cimp_cutoff, cimp_probes = MLH1_PROBES_450k, CIMP_CUTOFF_450k, CIMP_PROBES_450k
  elif array == '27k':
    mlh1_probes, cimp_cutoff, cimp_probes = MLH1_PROBES_27k, CIMP_CUTOFF_27k, CIMP_PROBES_27k
  else:
    logging.FATAL('Unrecognized array type %s', array)
    sys.exit(1)

  sys.stdout.write('Sample\tMLH1\tCIMP\t{}\t{}\n'.format('\t'.join(sorted(mlh1_probes)), '\t'.join(['\t'.join(['{}_{}'.format(gene, x) for x in sorted(cimp_probes[gene])]) for gene in sorted(cimp_probes)])))
  for f in open(files, 'r'):
    f = f.strip('\n')
    logging.info('processing %s...', f)
    mlh1_values = {}
    cimp_values = collections.defaultdict(list)
    probes = {}
    for line in open(f, 'r'):
      fields = line.strip('\n').split('\t')
      if fields[0] in mlh1_probes:
        mlh1_values[fields[0]] = fields[1]
      for gene in cimp_probes: # each cimp gene
        if fields[0] in cimp_probes[gene]: # relevant probe
          try:
            cimp_values[gene].append(float(fields[1]))
          except ValueError:
            logging.debug('skipping non-float value %s in probe %s', fields[1], fields[0])
          probes[fields[0]] = fields[1]

    sample = 'TCGA{}'.format(f.split('TCGA')[-1].replace('.gdc_hg38.txt', ''))
    sys.stdout.write('{}\t{:.2f}\t'.format(sample, sum([float(mlh1_values[x]) for x in mlh1_values]) / len(mlh1_values)))
    cimp_count = 0
    for gene in cimp_probes:
      if sum(cimp_values[gene]) / len(cimp_values[gene]) > cimp_cutoff:
        cimp_count += 1

    sys.stdout.write('{}\t'.format(cimp_count))

    # mlh1 probes
    sys.stdout.write('\t'.join(mlh1_values.get(x, 'NA') for x in sorted(mlh1_probes)))

    # cimp probes
    for gene in sorted(cimp_probes):
      for probe in sorted(cimp_probes[gene]):
        sys.stdout.write('\t{}'.format(probes.get(probe, 'NA')))
 
    sys.stdout.write('\n')

    logging.info('processing %s: done', f)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MLH1 methylation')
  parser.add_argument('--files', required=True, help='file containing filenames')
  parser.add_argument('--array', required=False, default='450k', help='450k or 27k')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.files, args.array)

