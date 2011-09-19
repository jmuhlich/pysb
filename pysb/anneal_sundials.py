import pysb.bng
import numpy 
import sympy 
import re 
import ctypes
import csv
import scipy.interpolate
from pysundials import cvode

# These set of functions set up the system for annealing runs
# and provide the runner function as input to annealing

def annlinit(model, reltol=1.0e-7, abstol=1.0e-11, nsteps = 1000, itermaxstep = None):
    '''
    must be run to set up the environment for annealing with pysundials
    '''
    # Generate equations
    pysb.bng.generate_equations(model)
    # Get the size of the ODE array
    odesize = len(model.odes)
    
    # init the arrays we need
    ydot = numpy.zeros(odesize) #dy/dt
    yzero = numpy.zeros(odesize)  #initial values for yzero
    
    # assign the initial conditions
    # FIXME: code outside of model shouldn't handle parameter_overrides 
    # FIXME: Species really should be a class with methods such as .name, .index, etc...
    for cplxptrn, ic_param in model.initial_conditions:
        override = model.parameter_overrides.get(ic_param.name)
        if override is not None:
            ic_param = override
        speci = model.get_species_index(cplxptrn)
        yzero[speci] = ic_param.value

    # initialize y with the yzero values
    y = cvode.NVector(yzero)
    numparams = len(model.parameters)
        
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

    #paramlist for annealing feeder function
    paramlist = []
    for i in range(0, numparams):
        # notice: p[i] ~ model.parameters[i].name ~ model.parameters[i].value
        data.p[i] = model.parameters[i].value
        paramlist.append(model.parameters[i].value)
    paramarray = numpy.asarray(paramlist)
    
    
    # if no sensitivity analysis is needed allocate the "p" array as a 
    # pointer array that can be called by sundials "f" as needed
    def f(t, y, ydot, f_data):
        data = ctypes.cast(f_data, PUserData).contents
        rhs_locals = {'y': y, 'p': data.p}
        for i in range(0,len(model.odes)):
            ydot[i] = eval(rhs_exprs[i], rhs_locals)
        return 0
    
    # CVODE STUFF
    # initialize the cvode memory object, use BDF and Newton for stiff
    cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
    # allocate the cvode memory as needed, pass the function and the init ys
    cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, reltol, abstol)
    # point the parameters to the correct array
    # if the params are changed later this does not need to be reassigned (???)
    cvode.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
    # link integrator with linear solver
    cvode.CVDense(cvode_mem, odesize)
    #stepsize
    if itermaxstep != None:
        cvode.CVodeSetMaxStep(cvode_mem, itermaxstep)

    #list of outputs
    xout = numpy.zeros(nsteps)
    yout = numpy.zeros([nsteps, odesize])

    #initialize the arrays
    #print "Initial parameter values:", y
    xout[0] = 0.0 #CHANGE IF NEEDED
    #first step in yout
    for i in range(0, odesize):
        yout[0][i] = y[i]
    
    return [f, rhs_exprs, y, ydot, odesize, data, xout, yout, nsteps, cvode_mem, yzero], paramarray


def annlodesolve(model, tfinal, envlist, params, useparams=None, tinit = 0.0, reltol=1.0e-7, abstol=1.0e-11):
    '''
    the ODE equation solver taylored to work with the annealing algorithm
    model: the model object
    envlist: the list returned from annlinit
    params: the list of parameters that are being optimized with annealing 
    useparams: the parameter number to which params[i] corresponds
    tinit: initial time
    reltol: relative tolerance
    abstol: absolute tolerance
    '''
    f = envlist[0]
    rhs_exprs = envlist[1]
    y = envlist[2]
    ydot = envlist[3]
    odesize = envlist[4]
    data = envlist[5]
    xout = envlist[6]
    yout = envlist[7]
    nsteps = envlist[8]
    cvode_mem = envlist[9]
    yzero = envlist[10]

    #set the initial values and params in each run
    #all parameters are used in annealing
    if useparams is None:
        for i in range(len(params)):
            data.p[i] = params[i]
    else:
        #only a subset of parameters are used for annealing
        for i in range(len(useparams)):
            #print "changing parameter", model.parameters[useparams[i]],"data.p", data.p[useparams[i]],"to", params[i]
            data.p[useparams[i]] = params[i]
        #for i, j in enumerate([x for x, y in enumerate(useparams) if y == 1]):
        #    data.p[j] = params[i]


    #reset initial concentrations
    y = cvode.NVector(yzero)

    # Reinitialize the memory allocations, DOES NOT REALLOCATE
    cvode.CVodeReInit(cvode_mem, f, 0.0, y, cvode.CV_SS, reltol, abstol)
    
    tadd = tfinal/nsteps

    t = cvode.realtype(tinit)
    tout = tinit + tadd
    
    #print "Beginning integration"
    #print "TINIT:", tinit, "TFINAL:", tfinal, "TADD:", tadd, "ODESIZE:", odesize
    print "Integrating Parameters:\n", params
    #print "y0:", yzero

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
    #print "Integration finished"

    #now deal with observables
    obs_names = [name for name, rp in model.observable_patterns]
    yobs = numpy.zeros([len(obs_names), nsteps])
    
    #sum up the correct entities
    for i, name in enumerate(obs_names):
        factors, species = zip(*model.observable_groups[name])
        yobs[i] = (yout[:, species] * factors).sum(1)

    #transpose yobs to make it easy to plot
    #yobs.T

    #merge the x and y arrays for easy analysis
    xyobs = numpy.vstack((xout, yobs))

    return (xyobs,xout,yout, yobs)

def read_csv_array(xpfname):
    """returns the first string and a numpy array from a csv set of data
    xpfname is a string with the file name"""

    #get the experimental data needed for annealing
    fp = open(xpfname, "r")

    reader = csv.reader(fp)
    templist = []
    #read in the lists
    for row in reader:
        templist.append(row)
    #remove empty spaces
    for i in range(0, len(templist)):
        templist[i] = filter(None, templist[i])
    headstring = templist.pop(0)
    #Now put these into a numpy array, assume they are all floats
    converters = tuple([float]*len(templist[0]))
    darray = numpy.zeros((len(templist), len(templist[0])))
    for i, item in enumerate(templist):
        darray[i] = numpy.asarray(templist[i], dtype=darray.dtype)    
   
    #transpose to ease analysis
    darray = darray.T
    return (darray, headstring)

def compare_data(xparray, simarray, xspairlist, xparrayvar=None):
    """Compares two arrays of different size and returns the X^2 between them.
    Uses the X axis as the unit to re-grid both arrays. 
    xparray: experimental data
    xparrayaxis: which axis of xparray to use for simulation
    simarray: simulation data
    simarrayaxis: which axis of simarray to use for simulation
    """
    # this expects arrays of the form array([time, measurement1, measurement2, ...])
    # the time is assumed to be roughly the same for both and the 
    # shortest time will be taken as reference to regrid the data
    # the regridding is done using a b-spline interpolation
    # xparrayvar shuold be the variances at every time point
    #
    # FIXME FIXME FIXME FIXME
    # This prob should figure out the overlap of the two arrays and 
    # get a spline of the overlap. For now just assume the simarray domain
    # is bigger than the xparray. FIXME FIXME FIXME
    #
    #rngmin = max(xparray[0].min(), simarray[0].min())
    #rngmax = min(xparray[0].max(), simarray[0].max())
    #rngmin = round(rngmin, -1)
    #rngmax = round(rngmax, -1)
    #print "Time overlap range:", rngmin,"to", rngmax
    
    ipsimarray = numpy.zeros(xparray.shape[1])
    objout = 0
   
    for i in range(len(xspairlist)):
        # create a b-spline of the sim data and fit it to desired range
        # import code
        # code.interact(local=locals())
        
        #some error checking
        #print "xspairlist length:", len(xspairlist[i])
        #print "xspairlist element type:", type(xspairlist[i])
        #print "xspairlist[i] elements:", xspairlist[i][0], xspairlist[i][1]
        assert type(xspairlist[i]) is tuple
        assert len(xspairlist[i]) == 2
        
        xparrayaxis = xspairlist[i][0]
        simarrayaxis = xspairlist[i][1]
        
        tck = scipy.interpolate.splrep(simarray[0], simarray[simarrayaxis])
        ipsimarray = scipy.interpolate.splev(xparray[0], tck) #xp x-coordinate values to extract from y splines
        
        # we now have x and y axis for the points in the model array
        # calculate the objective function
        #                        1
        # obj(t, params) = -------------(S_sim(t,params)-S_exp(t))^2
        #                  2*sigma_exp^2
        
        diffarray = ipsimarray - xparray[xparrayaxis]
        diffsqarray = diffarray * diffarray
        
        # assume a default .05 variance
        if xparrayvar is None:
            xparrayvar = numpy.ones(xparray.shape[1])
            xparrayvar = xparray[xparrayaxis]*.341 # 1 stdev w the experimental data... 
            xparrayvar = xparrayvar * xparrayvar

        xparrayvar = xparrayvar*2.0
        #numpy.seterr(divide='ignore')
        objarray = diffsqarray / xparrayvar

        # check for inf in objarray, they creep up when there are near zero or zero values in xparrayvar
        for i in range(len(objarray)):
            if numpy.isinf(objarray[i]) or numpy.isnan(objarray[i]):
                #print "CORRECTING NAN OR INF. IN ARRAY"
                # print objarray
                objarray[i] = 1e-100 #zero enough

        #import code
        #code.interact(local=locals())

        objout += objarray.sum()
        print "OBJOUT(%d,%d):%f  OBJOUT(CUM):%f"%(xparrayaxis, simarrayaxis, objarray.sum(), objout)

    print "OBJOUT(total):", objout
    return objout

def getgenparambounds(params, omag=3, N=1000.):
    # params must be a numpy array
    # from: http://projects.scipy.org/scipy/ticket/1126
    # The input-parameters "lower" and "upper" do not refer to global bounds of the
    # parameters space but to 'maximum' displacements in the MC scheme AND
    # in addition they determine the initial point!! 
    # The initial value that you provide for the parameter vector seems to have no influence.
    # This is how I call anneal with my desired functionality
    # p=[a,b,c] #my initial values
    # lb=array([a0,b0,c0]) #my lower bounds
    # ub=array([a1,b1,c1]) #my upper bounds
    # N=100 #determines the size of displacements
    # dx=(ub-lb)/N #displacements
    # lower=array(p)-dx/2 #the "lower bound" of the anneal routine
    # upper=array(p)+dx/2 #the "upper bound" of the anneal routine
    # f=lambda var: costfunction(p,lb,ub) #my cost function that is made very high if not lb < p < ub
    # pbest=scipy.optimize.anneal(f,p,lower=lower,upper=upper)
    # This ensures a MC search that starts of close to my initial value and makes steps of dx in its search.

    # set upper/lower bounds for generic problem
    ub = params * pow(10,omag)
    lb = params / pow(10,omag)
    
    # maximum displacements of the anneal routine
    dx = (ub-lb)/N
    lower = params - dx #/2
    upper = params + dx #/2
    
    return lb, ub, lower, upper


def annealfxn(params, useparams, time, model, envlist, xpdata, xspairlist, lb, ub, norm=False):
    ''' Feeder function for scipy.optimize.anneal
    '''
    # sample anneal call full model:
    # annlout = scipy.optimize.anneal(pysb.anneal_sundials.annealfxn, params, 
    #                                 args=(65000, model, envlist, xpdata, [(2,2), (3,3)], lb, ub), 
    #                                 lower=lower, upper=upper)
    # params: parameters to be optimized
    # lower,upper: arrays from get array function or something similar from getgenparambounds
    # lb, ub: lower bound and upper bound for function from getgenparambounds
    #
    # sample anneal call, optimization of some parameters
    #   annlout = scipy.optimize.anneal(pysb.anneal_sundials.annealfxn, smacprm, args=(smacnum, 25000, model, envlist, xpdata,
    #            [(2,2), (3,3)], lower=lower, upper=upper, full_output=1)
    #
    # sample anneal call, optimization for ALL parameters
    # 
    #

    if numpy.greater_equal(params, lb).all() and numpy.less_equal(params, ub).all():
        outlist = annlodesolve(model, time, envlist, params, useparams)
        # specify that this is normalized data
        if norm is True:
            print "Normalizing data"
            datamax = numpy.max(outlist[0], axis = 1)
            datamin = numpy.min(outlist[0], axis = 1)
            outlistnorm = ((outlist[0].T - datamin)/(datamax-datamin)).T
            # xpdata[0] should be time, get from original array
            outlistnorm[0] = outlist[0][0].copy()
            # xpdata here is normalized, and so is outlistnorm
            objout = compare_data(xpdata, outlistnorm, xspairlist)
        else:
            objout = compare_data(xpdata, outlist[0], xspairlist)
    else:
        print "======>VALUE OUT OF BOUNDS NOTED"
        temp = numpy.where((numpy.logical_and(numpy.greater_equal(params, lb), numpy.less_equal(params, ub)) * 1) == 0)
        for i in temp:
            print "======>",i, params[i]
        objout = 1.0e300 # the largest FP in python is 1.0e308, otherwise it is just Inf
    return objout

    





