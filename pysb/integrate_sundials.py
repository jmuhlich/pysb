import pysb.bng
import numpy, sympy, re, ctypes, sys
from pysundials import cvode, cvodes, nvecserial

#FIXME: add a odesolvesensedky function or add this function to odesenssolve

def odeinit(model):

    #Generate ODES from BNG
    # import code
    # code.interact(local=locals())
    # import pdb
    # pdb.set_trace()

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
    
    # Create the structure to hold the parameters when calling the function
    # This results in a generic "p" array
    class UserData(ctypes.Structure):
        _fields_ = [('p', cvode.realtype*numparams)] # parameters
    PUserData = ctypes.POINTER(UserData)
    data = UserData() 

    for i in range(0, numparams):
        # notice: p[i] ~ model.parameters[i].name ~ model.parameters[i].value
        data.p[i] = model.parameters[i].value

    # if no sensitivity analysis is needed allocate the "p" array as a 
    # pointer array that can be called by sundials "f" as needed
    def f(t, y, ydot, f_data):
        data = ctypes.cast(f_data, PUserData).contents
        rhs_locals = {'y': y, 'p': data.p}
        for i in range(0,len(model.odes)):
            ydot[i] = eval(rhs_exprs[i], rhs_locals)
        return 0

    return f, rhs_exprs, y, ydot, odesize, data

def odesolve(model, tfinal, nsteps = 100, tinit = 0.0, reltol=1.0e-7, abstol=1.0e-11):
    '''
    model: a model object from pysb
    tfinal: the final time to run the integration (in units of pysb model, usually second)
    _calspec: passed from calibrate_anneal to know which observable to use as calibration ref
    nsteps: the number of steps to take between time reports, default = 100
    tinit: the initial time-step, default = 0.0
    reltol: relative tolerance for sundials, default=1.0e-7
    abstol: absolute tolerance for sundials, default=1.0e-11
    '''
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
    xout = numpy.zeros(nsteps)
    yout = numpy.zeros([nsteps, odesize])

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

    #now deal with observables
    obs_names = [name for name, rp in model.observable_patterns]
    yobs = numpy.zeros([len(obs_names), nsteps])
    
    #sum up the correct entities
    for i, name in enumerate(obs_names):
        factors, species = zip(*model.observable_groups[name])
        yobs[i] = (yout[:, species] * factors).sum(1)

    #transpose yobs to make it easy to plot
    yout = yout.T
    return (xout,yobs,yout)

def odesenssolve(model, tfinal, nsteps = 100, tinit = 0.0, senslist=None,
                 sensmaglist=None, reltol=1.0e-7, abstol=1.0e-11):
    tadd = tfinal/nsteps

    SOMEFLAG = True
    if SOMEFLAG:
        f, rhs_exprs, y, ydot, odesize, data = odeinit(model)

    # initialize the cvode memory object
    cvodes_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)

    if senslist is None:
        #make a senslist for all parameters
        senslist = [n for n in range(0, len(model.parameters))]

    numsens = len(senslist)

    # the sensitivity function needs an array of the scaling factors for each parameter
    # for which sensitivity will be calculated
    # default a scale of "1" unless sensmaglist is passed
    if sensmaglist is None and senslist is None:
        sensmaglist = [1 for n in range(0, numsens)]
    elif sensmaglist is None and senslist is not None:
        #senslist was passed, assign mags of 1 to the items in senslist, 0 otherwise
        sensmaglist = [0 for n in range(0, numsens)]
        for n in senslist:
            sensmaglist[n] = 1
    else:
        print "something is really wrong with the SENSLIST or SENSMAGLIST"

    #FIXME: Should we remove the initial value sensitivities automatically?

    # set the sensitivity array
    yS = nvecserial.NVectorArray([([0]*odesize)]*numsens)

    # CVodeSensMalloc allocates and initializes memory for sensitivity computations
    cvodes.CVodeSensMalloc(cvodes_mem, numsens, cvodes.CV_STAGGERED, yS)

    # CVodeSetSensParams sets the parameters for the sensitivity function call
    print "SENSLIST:", senslist
    print "SENSMAGLIST:", sensmaglist

    cvodes.CVodeSetSensParams(cvodes_mem, data.p,
                              sensmaglist,
                              senslist)

    # point the user parameters to the correct array
    cvodes.CVodeSetFdata(cvodes_mem, ctypes.pointer(data))
    
    # allocate the cvodes memory as needed
    cvodes.CVodeMalloc(cvodes_mem, f, 0.0, y, cvodes.CV_SS, reltol, abstol)
    
    # link integrator with linear solver
    cvodes.CVDense(cvodes_mem, odesize)
    
    #list of outputs
    yout = numpy.zeros([nsteps, odesize])
    xout = numpy.zeros(nsteps)
    ysensout = numpy.zeros([odesize, nsteps, numsens])

    #initialize the arrays
    print "Initial parameter values:", y
    xout[0] = tinit
    for i in range(0, odesize):
        yout[0][i] = y[i]

    for n in range(0, numsens):
        for o in range(0, odesize):
            ysensout[i] = yS[n][o]

    t = cvode.realtype(tinit)
    tout = tinit + tadd

    print "Beginning integration, TINIT:", tinit, "TFINAL:", tfinal, "TADD:", tadd, "ODESIZE:", odesize
    for step in range(1, nsteps):
        ret = cvodes.CVode(cvodes_mem, tout, y, ctypes.byref(t), cvodes.CV_NORMAL)
        cvodes.CVodeGetSens(cvodes_mem, t, yS)
        
        if ret !=0:
            print "CVODES ERROR %i"%(ret)
            break

        xout[step]= tout
        for i in range(0, odesize):
            yout[step][i] = y[i]
            for j in range(0, numsens):
                #import code
                #code.interact(local=locals())
                ysensout[i][step][j] = yS[j][i] # yS[numsens][odesize]

            
        # increase the time counter
        tout += tadd
    print "Integration finished"

    #now deal with observables
    obs_names = [name for name, rp in model.observable_patterns]
    yobs = numpy.zeros([len(obs_names), nsteps])
    
    #sum up the correct entities
    for i, name in enumerate(obs_names):
        factors, species = zip(*model.observable_groups[name])
        yobs[0] = (yout[:, species] * factors).sum(1)

    #transpose yobs to make it easy to plot
    yobs.T
    return (xout, yobs, yout, ysensout)

# def read_csv_array(fp):
#     """returns the first string and a numpy array from a csv set of data"""
#     reader = csv.reader(fp)
#     templist = []
#     #read in the lists
#     for row in reader:
#         templist.append(row)
#     #remove empty spaces
#     for i in range(0, len(templist)):
#         templist[i] = filter(None, templist[i])
#     #Now put these into a numpy array
#     converters = (float, float, float, float)
#     headstring = templist.pop(0)
#     #assume all entries in the list are the same length
#     darray = numpy.zeros((len(templist), len(templist[0])))
#     for i, item in enumerate(templist):
#         darray[i] = numpy.asarray(templist[i], dtype=darray.dtype)    
#     return (headstring, darray)


# def compare_data(array0, array1, array0var=None, obspec=None):
#     """Compares two arrays of different size and returns the X^2 between them.
#     Uses the X axis as the unit to re-grid both arrays
#     obspec: passed by odeanneal to know which observable to use for comparison
#     """
#     # figure out the array shapes
#     # this expects arrays of the form array([time, measurements])
#     # the time is assumed to be roughly the same for both and the 
#     # shortest time will be taken as reference to regrid the data
#     # the regridding is done using a b-spline interpolation
#     # array0var shuold be the variances at every time point
#     #

#     # sanity checks
#     # make sure we are comparing the right shape arrays
#     arr0shape = array0.shape
#     arr1shape = array1.shape
    
#     if len(arr0shape) != len(arr1shape):
#         raise SystemExit("comparing arrays of different dimensions")
    
#     # get the time range where the arrays overlap
#     rngmin = max(array0[:,0].min(), array1[:,0].min)
#     rngmax = min(array0[:,0].max(), array1[:,0].max)
#     rngmin = round(rngmin, -1)
#     rngmax = round(rngmax, -1)
    
#     # use the experimental gridpoints from the reference array as
#     # the new gridset for the model array. notice the time range
#     # of the experiment has to be within the model range
#     #
#     iparray = numpy.zeros(array0.shape)
#     iparray[:,0] =  array0[:,0]
    
#     # now create a b-spline of the data and fit it to desired range
#     tck = scipy.interpolate.splrep(array1[:,0], array1[:,1])
#     iparray[:,1] = scipy.interpolate.splev(iparray[:,0], tck)
    
#     # we now have x and y axis for the points in the model array
#     # calculate the objective function
    
#     diffarray = array0[:,1] - iparray[:,1]
#     diffsqarray = diffarray * diffarray
    
#     # assume a default .05 variance
#     if array0var is None:
#         array0var = numpy.ones(iparray[:,1].shape)
#         array0var = array0var*.05
    
#     array0var = numpy.multiply(array0var,array0var)
#     array0var = array0var*0.5

#     objarray = diffsqarray * array0var
    
#     return objarray.sum()
    
# def odeanneal(refarray, params, model, tfinal, obspec=0):
#     '''
#     paramarray: array of parameters, from the model
#     eqns: dict of ode's
#     time: time to run the simulation
#     '''
#     # This is the feeder function for the anneal routine. It needs to:
#     # 1- run the integration for the model with parameters as vars
#     # 2- calculate the objective function and return it
#     # 3- return the difference
    
#     #FIXME: something about params here???? Say which params are being passed on by the fxn


#     # output contains xout, yobs, yout
#     output = odesolve(model, tfinal)
#     # specify which output is being compared!!
    
#     arrdiff = compare_data(refarray, output[obspec], )
    
#     return diff
