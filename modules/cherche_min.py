import sys
import numpy

file = open(sys.argv[1]);
values = file.readlines();
values = values[1:]
size = len(values)
cost_function = numpy.zeros(size)

for i in xrange(size):
	values[i] = values[i].split()
	cost_function[i] = float(values[-1])

file.close()
min = min(a)
indice = a.index(min)+1
print('Minimal value : {}. Iteration : {}.'.format(min,indice))
