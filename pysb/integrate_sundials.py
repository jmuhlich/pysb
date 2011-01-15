import pysb.bng
import numpy, sympy, re, ctypes
from pysundials import cvode

#FIXME: add a odesolvesensedky function or add this function to odesenssolve

def odeinit(model, senslist=None):

    #Generate ODES from BNG
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
    numparams = len(model.parameters)
        
    #print "initial parameter values:\n", y

    # make a dict of ydot functions. notice the functions are in this namespace.
    # replace the kxxxx constants with elements from the params array
    rhs_exprs = []
    for i in range(0,odesize):
        # first get the function string from sympy, replace the the "sN" with y[N]
        tempstring = re.sub(r's(\d+)', lambda m: 'y[%s]'%(int(m.group(1))), str(model.odes[i]))
        # now replace the constants with 'p' array names; cycle through the whole list
        for j in range(0, numparams):
            tempstring = re.sub('(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_])'%(model.parameters[j].name),
                                'p[%d]'%(j), tempstring)

        # make a list of compiled rhs expressions which will be run by the integrator
        # use the ydots to build the function for analysis
        # (second arg is the "filename", useful for exception/debug output)
        rhs_exprs.append(compile(tempstring, '<ydot[%s]>' % i, 'eval'))

    # senslist defined only if sensitivity list passed by function calling ODEinit. 
    # if senslist empty then the allocate space for sensitivity of all parameters.
    if senslist:
        #assign the ctypes object for sensitivity analysis
        if len(senslist) is 0:
            sensnum = len(model.parameters)
            senslist = [n for n in range(0, sensnum)]
        else:
            sensnum = len(senslist)

        #create the structure to hold the parameters when calling the function
        # This results in a generic "p" array
        class UserData(ctypes.Structure):
            _fields_ = [('p', cvode.realtype*sensnum)] # parameters
        PUserData = ctypes.POINTER(UserData)
        data = UserData() 

        # Store the parameter values in the cvodes array
        for i in range(0, sensnum):
            # notice: p[i] ~ model.parameters[i].name ~ model.parameters[i].value
            data.p[i] = model.parameters[senslist[i]].value

        #Finally define the function that will be integrated
        def f(t, y, ydot, f_data):
            data = ctypes.cast(f_data, PUserData).contents
            for i in range(0,len(model.odes)):
                ydot[i] = eval(rhs_exprs[i], None, {'y': y, 'p': data.p})
            return 0
    else:
        #create the structure to hold the parameters when calling the function
        # This results in a generic "p" array
        class UserData(ctypes.Structure):
            _fields_ = [('p', cvode.realtype*numparams)] # parameters
        PUserData = ctypes.POINTER(UserData)
        data = UserData() 

        # if no sensitivity analysis is needed allocate the "p" array as a 
        # pointer array that can be called by sundials "f" as needed
        for i in range(0, numparams):
            # notice: p[i] ~ model.parameters[i].name ~ model.parameters[i].value
            data.p[i] = model.parameters[i].value
        def f(t, y, ydot, f_data):
            data = ctypes.cast(f_data, PUserData).contents
            rhs_locals = {'y': y, 'p': data.p}
            for i in range(0,len(model.odes)):
                ydot[i] = eval(rhs_exprs[i], rhs_locals)
            return 0

    return f, rhs_exprs, y, ydot, odesize, data

def odesolve(model, tfinal, nsteps = 100, tinit = 0.0, reltol=1.0e-8, abstol=1.0e-12):
    tadd = tfinal/nsteps

    SOMEFLAG = True
    if SOMEFLAG:
        f, rhs_exprs, y, ydot, odesize, data = odeinit(model)

    # initialize the cvode memory object
    cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
    
    # allocate the cvode memory as needed
    cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, reltol, abstol)
    
    # point the parameters to the correct array
    cvode.CVodeSetFdata(cvode_mem, ctypes.pointer(data))

    # link integrator with linear solver
    cvode.CVDense(cvode_mem, odesize)
    
    #list of outputs
    yout = numpy.zeros([nsteps, odesize])
    xout = numpy.zeros(nsteps)

    #initialize the arrays
    print "Initial parameter values:", y
    xout[0] = tinit
    for i in range(0, odesize):
        yout[0][i] = y[i]

    t = cvode.realtype(tinit)
    tout = tinit + tadd

    print "Beginning integration, TINIT:", tinit, "TFINAL:", tfinal, "TADD:", tadd, "ODESIZE:", odesize
    for step in range(1, nsteps):

        ret = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
        
        if ret !=0:
            print "CVODE ERROR %i"%(ret)
            break

        xout[step]= tout
        for i in range(0, odesize):
            yout[step][i] = y[i]

        # increase the time counter
        tout += tadd
    print "Integration finished"

    return (xout,yout)

def odesenssolve(model, tfinal, senslist=None, reltol=1.0e-8, abstol=1.0e-12):
    SOMEFLAG = True
    if SOMEFLAG:
        f, rhs_exprs, y, ydot, odesize, data = odeinit(model, senslist)
    
    # initialize the cvode memory object
    cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
    
    # allocate the cvode memory as needed
    cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, reltol, abstol)
    
    # link integrator with linear solver
    cvode.CVDense(cvode_mem, odesize)
    
    #list of youts
    yout = []
    for i in range(0, odesize+1): #leave one for the timestamp
        yout.append([])

    t = cvode.realtype(0.0)
    iout = 0 #initial time
    tout = 0.1 #time of next integration
    
    print "Beginning integration"
    while iout < tfinal:
        ret = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
        
        if ret !=0:
            print "CVODE ERROR %i"%(ret)
            break
        yout[0].append(tout)
        for i in range(1, odesize):
            yout[i].append(y[i])

        # increase the while counter
        iout += 1
    print "Integration finished"

    return yout
