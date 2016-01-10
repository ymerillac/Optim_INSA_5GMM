import sys
import numpy

file = open(sys.argv[1])
values = file.readlines()
values = values[1:]
size = len(values)
cost_function = numpy.zeros(size)

for i in xrange(size):
	tmp = values[i].split()
	if float(tmp[-2]) < 3e-3 :
		cost_function[i] = 1
	else :
		cost_function[i] = float(tmp[-2])
	constraint_function = float(tmp[-1])
	if constraint_function > 0:
		cost_function[i] = max(cost_function)+0.1

file.close()
min = min(cost_function)
indice = cost_function.argmin()+1
print('Minimal value : {}. Iteration : {}.'.format(min,indice))
