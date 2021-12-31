#!/usr/bin/env python

import csv
import sys

#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
#0       NM_001376549.1  chr1    +       66999043        67216822        67000041        67208778        23      66999043,66999928,67091529,67098752,67101626,67105459,67108492,67109226,67136677,67137626,67138963,67142686,67145360,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755,       66999090,67000051,67091593,67098777,67101698,67105516,67108547,67109402,67136702,67137678,67139049,67142779,67145435,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,67216822, 0    SGIP1    cmpl    cmpl    -1,0,1,2,0,0,0,1,0,1,2,1,1,1,0,1,1,2,2,0,2,1,1,
def make_canonical(ifh, ofh):
  best = {}
  for r in csv.DictReader(ifh, delimiter='\t'):
    tx = r['name']
    txlen = int(r['txEnd']) - int(r['txStart'])
    gene = r['name2']
    if gene not in best:
      best[gene] = (tx, txlen)
    elif txlen > best[gene][1]:
      best[gene] = (tx, txlen)

  # results
  o = csv.DictWriter(ofh, delimiter='\t', fieldnames=['gene', 'longest', 'length'])
  o.writeheader()
  for b in sorted(best):
    o.writerow({'gene': b, 'longest': best[b][0], 'length': best[b][1]})

if __name__ == '__main__':
  make_canonical(sys.stdin, sys.stdout)
