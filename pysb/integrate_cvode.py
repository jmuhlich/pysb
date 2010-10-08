import pysb.bng
import numpy
from pysundials import cvode
import sympy
import re


def odesolve(model, t):
    # Generate the ODES (Jah sucks at comments)
    pysb.bng.generate_equations(model)

    # Get the size of the ODE array
    odesize = len(model.odes)
    
    #make a function that can be passed to cvode
    #init the arrays we need
    y = numpy.zeros(len(model.odes))
    ydot = y.copy() 
    y0 = y.copy()
    
    # assign the initial conditions
    # FIXME code outside of model shouldn't have to handle parameter_overrides (same for initial_conditions below)
    # CFL: Species really should be a class with methods such as .name, .index, etc
    for cplxptrn, ic_parm in model.initial_conditions:
        override = model.parameter_overrides.get(ic_parm.name)
        if override is not None:
            ic_parm = override
        speci = model.get_species_index(cplxptrn)
        y0[speci] = ic_parm.value
    
    # define function to integrate and assign RHS from BNG
    def f(t, y, ydot, f_data):
        for i in 
        exec "%s = %f" % (model.parameters[0].name, model.parameters[0].value)



    # Get the constants from the model
    c_code_consts = '\n'.join(['float %s = %e;' % (p.name, p.value) for p in model.parameters])


    c_code_eqs = '\n\t'.join(['ydot[%d] = %s;' % (i, sympy.ccode(model.odes[i])) for i in range(len(model.odes))])
    c_code_eqs = re.sub(r's(\d+)', lambda m: 'y[%s]' % (int(m.group(1))), c_code_eqs)
    c_code = c_code_consts + '\n\n' + c_code_eqs

    y0 = numpy.zeros((len(model.odes),))
        
    def rhs(y, t):
        ydot = y.copy()  # seems to be the fastest way to get an array of the same size?
        inline(c_code, ['y', 'ydot']); # sets ydot as a side effect
        return ydot

    nspecies = len(model.species)
    obs_names = [name for name, rp in model.observable_patterns]
    rec_names = ['__s%d' % i for i in range(nspecies)] + obs_names
    yout = numpy.ndarray((len(t), len(rec_names)))

    # perform the actual integration
    yout[:, :nspecies] = odeint(rhs, y0, t)

    for i, name in enumerate(obs_names):
        factors, species = zip(*model.observable_groups[name])
        yout[:, nspecies + i] = (yout[:, species] * factors).sum(1)

    dtype = zip(rec_names, (yout.dtype,) * len(rec_names))
    yrec = numpy.recarray((yout.shape[0],), dtype=dtype, buf=yout)
    return yrec
