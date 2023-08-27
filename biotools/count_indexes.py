#!/usr/bin/env python

import collections
import gzip
import sys

i = collections.defaultdict(int)
for l in sys.stdin:
  if l.startswith('@'):
    ix = l.strip()[-8:]
    if ix not in i:
      sys.stderr.write('new index: {}\n'.format(ix))
    i[ix] += 1

for ix in sorted(i.keys()):
  sys.stdout.write('{}\t{}\n'.format(ix, i[ix]))
