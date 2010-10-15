import pysb.bng
import numpy, sympy, re
from pysundials import cvode

def odesolve(model, output):
    # Generate the ODES (Jah sucks at comments)
    pysb.bng.generate_equations(model)

    # Get the size of the ODE array
    odesize = len(model.odes)
    
    # init the arrays we need
    ydot = numpy.zeros(odesize) #dy/dt
    y0 = y.copy() # initial values for y (bound)

    # assign the initial conditions
    # FIXME: code outside of model shouldn't handle parameter_overrides 
    # FIXME: Species really should be a class with methods such as .name, .index, etc... jah is good at this
    for cplxptrn, ic_parm in model.initial_conditions:
        override = model.parameter_overrides.get(ic_parm.name)
        if override is not None:
            ic_parm = override
        speci = model.get_species_index(cplxptrn)
        y0[speci] = ic_parm.value
    
    # get parameters from BNG
    # get a key:value parameters dict. notice these are local values
    for i in range(0, len(model.parameters)):
        key = model.parameters[i].name
        if "_0" in key:
            continue
        else:
            exec "%s = %f" % (model.parameters[i].name, model.parameters[i].value)

    # make a dict of ydot functions. notice the functions are in this namespace.
    funcs = {}
    for i in range(0,len(model.odes)):
        exec "def _ydot%d(y): return %s" % (i, re.sub(r's(\d+)', lambda m: 'y[%s]'%(int(m.group(1))), str(model.odes[i])))
        funcs[i] = eval("_ydot%d"%(i))
    
    # use the _ydots to build the function for analysis
    # FIXME: the best way i could think of not doing an "exec" or an "eval" in each loop was to
    #        map a dict to a function. although I love python I miss C pointers. Perhaps ctypes?
    def f(t, y, ydot, f_data):
        for i in range(0,len(model.odes)):
            ydot[i] = funcs[i](y)
        return 0
    
    # initialize y
    y = cvode.Nvector(y0)

    # initialize the cvode memory object
    cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
    
    # allocate the cvode memory as needed
    cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, 1.0e-8, 1.0e-12)
    
    # link integrator with linear solver
    cvode.CVDense(cvode_mem, odesize)
    
    #list of outputs
    output = []
    for i in range(0, odesize+1): #leave one for the timestamp
        output.append([])

    iout = 0 #initial time
    tout = 0.1 #time of next integration
    
    print "Beginning integration"
    while iout < 10000:
        ret = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
        
        if ret !=0:
            print "CVODE ERROR %i"%(ret)
            break
        output[0].append(tout)
        for i in range(1, odesize):
            output[i].append(y[i])

        # increase the while counter
        iout += 1
    print "Integration finished"

    return output
