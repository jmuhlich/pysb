import pysb
import pysb.bng
from pylab import *
from scipy.integrate import odeint
from scipy.weave import inline
import sympy
import re
import time
from pprint import pprint

from earm_1_0 import model


pysb.bng.generate_equations(model)

#print "EQUATIONS:"
#for i, expr in enumerate(model.odes):
#    print 'ds%d/dt = %s' % (i + 1, expr)
#print

y_syms = [ sympy.Symbol('s%d' % n) for n in range(1,len(model.odes)+1) ]
param_subs = dict([ (sympy.Symbol(p.name), p.value) for p in model.parameters ])

c_code_consts = '\n'.join(['float %s = %e;' % (p.name, p.value) for p in model.parameters])
c_code_eqs = '\n'.join(['ydot[%d] = %s;' % (i, sympy.ccode(model.odes[i])) for i in range(len(model.odes))])
c_code_eqs = re.sub(r's(\d+)', lambda m: 'y[%s]' % (int(m.group(1))-1), c_code_eqs)
c_code = c_code_consts + '\n\n' + c_code_eqs

def rhs(y, t):
    global rhs_time

    tic = time.time();
    ydot = y.copy()  # seems to be the fastest way to get an array of the same size?
    inline(c_code, ['y', 'ydot']); # sets ydot as a side effect
    toc = time.time();
    rhs_time += toc - tic;
    return ydot

y0_explicit = [model.parameter(m.name + '_0').value for m in model.monomers if model.parameter(m.name + '_0') is not None]
y0 = concatenate((y0_explicit, zeros(len(model.odes)-len(y0_explicit))))

t = linspace(0, 6*3600, 1000)



rhs_time = 0
t_start = time.time()

yout = odeint(rhs, y0, t)

t_end = time.time()

print "total time:", t_end - t_start
print "RHS CPU time:", rhs_time
print "fraction of total time in RHS:", rhs_time/(t_end - t_start)


# reproduce Fig. 4B (right side) from plos biol paper

yplot = zeros((len(t), 3))
#yplot[:,0] = yout[:,(28-1,32-1,33-1)].sum(1)
#yplot[:,1] = yout[:,36-1]
#yplot[:,2] = yout[:,(52-1,54-1)].sum(1)
yplot = yout[:,(9,8,15)]
yplot /= yplot.max(0)
yplot = 1-yplot

tp = t / 3600
p = plot(tp, yplot[:,0], 'b', tp, yplot[:,1], 'y', tp, yplot[:,2], 'r')
figlegend(p, ('tBid', 'CPARP', 'cSmac'), 'upper left')
a = gca()
a.set_ylim((-.05, 1.05))

#show()
