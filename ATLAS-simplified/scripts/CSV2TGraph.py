#!/usr/bin/env python

import sys
from array import array
from math import fabs

try:
  input_file=sys.argv[1]
except:
  print 'Error: no input file provided...'; quit()

x=array('f');
y=array('f');
_colunn = 3; # specify column number to read from

print 'Read data from csv file'
for line in [line.strip() for line in open(str(input_file))]:
  if '#' in line: continue
  if len(line)==0: continue
  data = line.split(',')
  x.append(float(data[0]))
  y.append(float(data[_colunn]))
  print 'add (',data[0],',',data[_colunn],')'

n=len(x)
print('if (x<%2.4f) return %2.6f;'%(x[0],y[0]))
for i in range(1,n):
  m = (y[i] - y[i-1])/(x[i] - x[i-1])
  print('else if (x<%2.4f) return %2.6f%s%2.8f*x;'%(x[i],y[i-1] - m*x[i-1],(' + ' if (m>0) else ' - '),fabs(m)))
print('else return %2.6f;'%y[n-1])
