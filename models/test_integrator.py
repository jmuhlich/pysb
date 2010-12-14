import pysb
import pysb.bng
from pylab import *
from scipy.integrate import odeint
import sympy
import time

pysb.bng.pkg_path = '/usr/local/share/bionetgen'

from earm_1_0 import model, bng_content

ydot = pysb.bng.generate_equations(bng_content)

#print "EQUATIONS:"
#for s in sorted(ydot.keys()):
#    print 'ds%d/dt = %s' % (s, ydot[s])
#print

y_syms = [ sympy.Symbol('s%d' % n) for n in range(1,len(ydot)+1) ]
param_subs = dict([ (sympy.Symbol(p.name), p.value) for p in model.parameters ])
ydot_list = [ ydot[k] for k in sorted(ydot.keys()) ]

def rhs(y, t0):
    global rhs_time

    tic = time.time();
    subs = param_subs.copy()
    subs.update([ (y_syms[i], y[i]) for i in range(len(y)) ])
    ret = [ exp.evalf(subs=subs) for exp in ydot_list ]
    toc = time.time();
    rhs_time += toc - tic;
    return ret

y0_explicit = [model.parameter(m.name + '_0').value for m in model.monomers if model.parameter(m.name + '_0') is not None]
y0 = concatenate((y0_explicit, zeros(len(ydot)-len(y0_explicit))))

#t = linspace(0, 6*3600, 6*60+1)
t = linspace(0, 400, 401)  # python sympy rhs is way too slow to simulate 6 hours



rhs_time = 0
t_start = time.time()

yout = odeint(rhs, y0, t)

t_end = time.time()



print "total time:", t_end - t_start
print "RHS CPU time:", rhs_time
print "fraction of total time in RHS:", rhs_time/(t_end - t_start)

p = plot(t, yout)
l = figlegend(p, y_syms, 'upper right')
a = gca()
a.set_yscale('log')
a.set_ylim(1e-80,1e6)

show()
