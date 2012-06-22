#!/usr/bin/python

from numpy import arange, cos, exp
import pylab

x = arange(-100.0, 100.0, 0.01)
y = exp(x)

pylab.plot(x, y)
pylab.show()