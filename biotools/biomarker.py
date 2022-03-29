#!/usr/bin/env python
'''
  chi-square on common variants
'''

import argparse
import collections
import csv
import logging
import sys

import cyvcf2
import statx
import statx.chi

SAMPLE_COL='Filename'
GROUP_COL='pks_positive'
PREFIX="vcfs/"
SUFFIX=".intersect.pass.filter.vcf.gz"

def main(groups):
  logging.info('reading stdin...')
  variants = collections.defaultdict(dict)
  hotspots = set()
  groups = collections.defaultdict(int)
  genes = {}
  for row in csv.DictReader(sys.stdin, delimiter='\t'):
    s = row[SAMPLE_COL]
    g = row[GROUP_COL]
    groups[g] += 1

    logging.info('processing %s...', s)
    
    vcf_in = cyvcf2.VCF('{}{}{}'.format(PREFIX, s, SUFFIX))
    for line, variant in enumerate(vcf_in):
      key = '{}/{}/{}/{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
      if g not in variants[key]:
        variants[key][g] = 0
      variants[key][g] += 1

      gene_pos = '{}/{}'.format(variant.CHROM, variant.POS)
      if gene_pos not in genes:
        #Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL
        gene = variant.INFO['CSQ'].split(',')[0].split('|')[5]
        genes[gene_pos] = gene
        logging.debug('%s -> %s', gene_pos, gene)

      if sum([variants[key][x] for x in variants[key]]) > 1:
        logging.debug('%s is a hotspot', key)
        hotspots.add(key)

    logging.info('processing %s: %i variants added', s, line)
    logging.info('%i hotspots so far', len(hotspots))
    

  logging.info('%i hotspots from %s', len(hotspots), groups)
  sys.stdout.write('Variant\tGene\tpvalue\t{}\n'.format('\t'.join(sorted(list(groups)))))
  for h in list(hotspots):
    # generate a p-value 'SPS-HMG': 7, 'SPS-IMG': 5, 'CRC-HMG': 5, 'CRC-IMG': 7, 'CRC-NL': 2
    #target = ('SPS-HMG',)
    tp = sum([variants[h].get(x, 0) for x in groups]) # target with condition
    fn = sum([groups[x] for x in groups]) - tp # target without condition
    fp = sum([variants[h].get(x, 0) for x in groups if x not in groups]) # non-target with condition
    tn = sum([groups[x] for x in groups if x not in groups]) - fp # non-target without condition
    result = statx.chi.fisher([tp, fp, fn, tn])

    gene = genes['/'.join(h.split('/')[:2])]
    sys.stdout.write('{}\t{}\t{:.6f}\t{}\n'.format(h, gene, result['pvalue'], '\t'.join([str(variants[h].get(g, 0)) for g in sorted(list(groups))])))
  
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Meh')
  parser.add_argument('--groups', nargs='+', required=True, help='positive groups')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.groups)

