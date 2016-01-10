import sys
import numpy

file = open(sys.argv[1])
values = file.readlines()
values = values[1:]
size = len(values)
cost_function = numpy.zeros(size)

for i in xrange(size):
	tmp = values[i].split()
	cost_function[i] = float(tmp[-1])

file.close()
min = min(cost_function)
indice = cost_function.argmin()+1
print('Minimal value : {}. Iteration : {}.'.format(min,indice))
