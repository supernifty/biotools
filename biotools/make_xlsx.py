#!/usr/bin/env python
# target, tsvs
# usage: make_xlsx.py --segment 1 out.cidr/all.xlsx out.cidr/all.*.blah.tsv
# usage: make_xlsx.py --segment 0 out.cidr/all.xlsx out.cidr/*.blah.tsv
# usage: make_xlsx.py --segment 2 out.cidr/all.xlsx out.cidr/all.foos.*.blah.tsv

import argparse
import logging
import sys

import pandas as pd

def make_xlsx(target, tsvs, delimiter, segment=2):
  logging.info('writing to %s...', target)
  writer = pd.ExcelWriter(target, engine='xlsxwriter')
  for tsv in tsvs:
    logging.info('adding %s...', tsv)
    df = pd.read_csv(tsv, sep=delimiter)
    #sheet_name='-'.join(tsv.split('/')[-1].split('.')[:2])[:31] # first 31 chars of filename
    sheet_name=tsv.split('/')[-1].split('.')[segment][:31] # first 31 chars of filename
    df.to_excel(writer, index=False, sheet_name=sheet_name)

  #writer.save()
  writer._save()
  writer.close()

  sys.stderr.write('done\n')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--segment', required=False, type=int, default=0, help='which part of filename to use for excel tab name')
  parser.add_argument('--target', required=True, help='write to this file')
  parser.add_argument('--delimiter', required=False, default=',', help='input file delimiter')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('csvs', nargs='+', help='csv files to merge')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  make_xlsx(args.target, args.csvs, args.delimiter, args.segment)
