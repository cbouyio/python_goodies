#!/usr/bin/env python

import sys

n = int(sys.argv[1])
name = sys.argv[2]
picName = sys.argv[3]

for i in range(n) :
  i = i + 1
  dir = "%s%02i/" % (name, i,)
  print '%i) Repetition %i\n```{r %s%i, rgl=TRUE, fig.width=15, fig.hight=15}\nplotFiberInteractions("%s")\n```\n' % (i, i, picName, i, dir,)

