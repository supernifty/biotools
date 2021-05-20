#!/usr/bin/env python
# target, tsvs
# usage: make_xlsx.py --segment 1 out.cidr/all.xlsx out.cidr/all.*.blah.tsv
# usage: make_xlsx.py --segment 0 out.cidr/all.xlsx out.cidr/*.blah.tsv
# usage: make_xlsx.py --segment 2 out.cidr/all.xlsx out.cidr/all.foos.*.blah.tsv

import sys

import pandas as pd

def make_xlsx(target, tsvs, segment=2):
  sys.stderr.write('writing to {}...\n'.format(target))
  writer = pd.ExcelWriter(target, engine='xlsxwriter')
  for tsv in tsvs:
    sys.stderr.write('adding {}...\n'.format(tsv))
    df = pd.read_csv(tsv, sep='\t')
    #sheet_name='-'.join(tsv.split('/')[-1].split('.')[:2])[:31] # first 31 chars of filename
    sheet_name=tsv.split('/')[-1].split('.')[segment] # first 31 chars of filename
    df.to_excel(writer, index=False, sheet_name=sheet_name)

  writer.save()
  writer.close()

  sys.stderr.write('done\n')

if __name__ == '__main__':
  if sys.argv[1] == '--segment':
    make_xlsx(sys.argv[3], sys.argv[4:], int(sys.argv[2]))
  else:
    make_xlsx(sys.argv[1], sys.argv[2:])
