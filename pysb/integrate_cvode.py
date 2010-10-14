import pysb.bng
import numpy
from pysundials import cvode
import sympy
import re



def odeinit(model):
    # Generate the ODES (Jah sucks at comments)
    pysb.bng.generate_equations(model)

    # Get the size of the ODE array
    odesize = len(model.odes)
    
    #make a function that can be passed to cvode
    #init the arrays we need
    y = numpy.zeros(len(model.odes)) #changing values for integration y[0]...y[n]
    ydot = y.copy() # dy/dt
    y0 = y.copy() # initial values for y (bound)
    
    # assign the initial conditions
    # FIXME code outside of model shouldn't handle parameter_overrides 
    # CFL: Species really should be a class with methods such as .name, .index, etc
    for cplxptrn, ic_parm in model.initial_conditions:
        override = model.parameter_overrides.get(ic_parm.name)
        if override is not None:
            ic_parm = override
        speci = model.get_species_index(cplxptrn)
        y0[speci] = ic_parm.value
    
    # get parameters from BNG
    keyvals = {}
    # first get a key:value parameters dict
    for i in range(0, len(model.parameters)):
        key = model.parameters[i].name
        if "_0" in key:
            continue
        else:
            exec "%s = %f" % (model.parameters[i].name, model.parameters[i].value)

    # make a set of ydot functions. notice the functions are in this namespace.
    for i in range(0,len(model.odes)):
        exec "def _ydot%s(y): return %s" % (i, str(re.sub(r's(\d+)', lambda m: 'y[%s]' % (int(m.group(1))), model.odes[i])))

    #use the _ydots to build the function for analysis
    #FIXME: this is probably going to be RIDICULOUSLY slow to do an eval on each integration step... 
    #I can't think of a better way to do this unless we do this into a class and somehow
    #declare the functions at a lower namespace... OR declare functions right from time we read the BNG code...
    def f(t, y, ydot, f_data):
        for i in range(0,len(model.odes)):
            ydot[i] = eval("_ydot%d(y)"%(i))
        return 0
    


    # Get the constants from the model

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
