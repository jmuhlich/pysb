import pysb.bng
import numpy, sympy, re, ctypes
from pysundials import cvode

#FIXME: this should have an ODEINIT section function and a few functions, an
#       odesolve, an odesolveSens, and odesolvesensedky or something like this

def odeinit(model, senslist=None):

    # senslist is defined only if the sensitivity list is passed
    # by the function calling ODEinit. This is just a sanity check.
    # if the sensitivity list is empty then allocate space for 
    # sensitivity analysis of the whole system.
    if senslist:
        #assign the ctypes object for sensitivity analysis
        if len(senslist) is 0:
            sensnum = len(model.parameters)
            senslist = [n for n in range(0, sensnum)]
        else:
            sensnum = len(senslist)
        
        class UserData(ctypes.Structure):
            _fields_ = [('p', cvodes.realtype*sensnum)] # parameters
        PUserData = ctypes.POINTER(UserData)
        data = UserData() 

        for i in range(0, sensnum):
            data.p[i] = model.parameters[senslist[i]].value

    #Generate ODES
    pysb.bng.generate_equations(model)
    
    # Get the size of the ODE array
    odesize = len(model.odes)
    
    # init the arrays we need
    ydot = numpy.zeros(odesize) #dy/dt
    yzero = numpy.zeros(odesize)  #initial values for yzero
    
    # assign the initial conditions
    # FIXME: code outside of model shouldn't handle parameter_overrides 
    # FIXME: Species really should be a class with methods such as .name, .index, etc...
    for cplxptrn, ic_parm in model.initial_conditions:
        override = model.parameter_overrides.get(ic_parm.name)
        if override is not None:
            ic_parm = override
        speci = model.get_species_index(cplxptrn)
        yzero[speci] = ic_parm.value

    # initialize y with the yzero values
    y = cvode.NVector(yzero)
    
    # get parameters from BNG
    # get a key:value parameters dict. notice these are local values
    numparams = len(model.parameters)
    params = numpy.zeros(numparams)
    for i in range(0, numparams):
        # notice: params[i] ~ model.parameters[i].name ~ model.parameters[i].value
        # IN THAT ORDER
        params[i] = model.parameters[i].value
        
    print "initial parameter values:\n", y

    # make a dict of ydot functions. notice the functions are in this namespace.
    # replace the kxxxx constants with elements from the params array
    funcs = {}
    for i in range(0,odesize):
        tempstring = re.sub(r's(\d+)', lambda m: 'y[%s]'%(int(m.group(1))), str(model.odes[i]))
        for j in range(0, numparams):
            tempstring = tempstring.replace(model.parameters[j].name, "p[%d]"%j)
        def tempfunc(y): return eval(tempstring)
        funcs[i] = tempfunc

    # use the ydots to build the function for analysis
    # FIXME: the best way i could think of not doing an "exec" or an "eval" in each loop was to
    #        map a dict to a function. although I love python I miss C pointers. Perhaps ctypes?
    def f(t, y, ydot, f_data):
        for i in range(0,len(model.odes)):
            ydot[i] = funcs[i](y)
        return 0
    
    return f, funcs, y, ydot, odesize



def odesolve(model, tfinal):
    SOMEFLAG = True
    if SOMEFLAG:
        f, funcs, y, ydot, odesize = odeinit(model)
    
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

    t = cvode.realtype(0.0)
    iout = 0 #initial time
    tout = 0.1 #time of next integration
    
    print "Beginning integration"
    while iout < tfinal:
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

def odesenssolve(model, tfinal, senslist=None, reltol=1.0e-8, abstol=1.0e-12):
    SOMEFLAG = True
    if SOMEFLAG:
        f, funcs, y, ydot, odesize = odeinit(model, senslist)
    
    # initialize the cvode memory object
    cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
    
    # allocate the cvode memory as needed
    cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, reltol, abstol)
    
    # link integrator with linear solver
    cvode.CVDense(cvode_mem, odesize)
    
    #list of outputs
    output = []
    for i in range(0, odesize+1): #leave one for the timestamp
        output.append([])

    t = cvode.realtype(0.0)
    iout = 0 #initial time
    tout = 0.1 #time of next integration
    
    print "Beginning integration"
    while iout < tfinal:
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
