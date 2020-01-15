#!/usr/bin/env python

import sys

PROBES=set(['cg23658326', 'cg11600697', 'cg21490561'])

values = {}
for line in sys.stdin:
  fields = line.strip('\n').split('\t')
  if fields[0] in PROBES:
    values[fields[0]] = fields[1]

sys.stdout.write('{:.2f}\t'.format(sum([float(values[x]) for x in values]) / len(values)))
sys.stdout.write('\t'.join(values.get(x, 'NA') for x in PROBES))
sys.stdout.write('\n')
