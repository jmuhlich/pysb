import pysb.bng
import numpy, sympy, re, ctypes
import csv
from pysundials import cvode

# These set of functions set up the system for annealing runs
# and provide the runner function as input to annealing

def annealinit(model, reltol=1.0e-7, abstol=1.0e-11, nsteps = 100):
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
    for cplxptrn, ic_parm in model.initial_conditions:
        override = model.parameter_overrides.get(ic_parm.name)
        if override is not None:
            ic_parm = override
        speci = model.get_species_index(cplxptrn)
        yzero[speci] = ic_parm.value

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
    
    # CVODE STUFF
    # initialize the cvode memory object, use BDF and Newton for stiff
    cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
    # allocate the cvode memory as needed, pass the function and the init ys
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
    xout[0] = 0.0 #CHANGE IF NEEDED
    #first step in yout
    for i in range(0, odesize):
        yout[0][i] = y[i]

    return [f, rhs_exprs, y, ydot, odesize, data, xout, yout, nsteps, cvode_mem]


def annealodesolve(model, tfinal, initlist, tinit = 0.0):
    f = initlist[0]
    rhs_exprs = initlist[1]
    y = initlist[2]
    ydot = initlist[3]
    odesize = initlist[4]
    data = initlist[5]
    xout = initlist[6]
    yout = initlist[7]
    nsteps = initlist[8]
    cvode_mem = initlist[9]

    tadd = tfinal/nsteps

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
        yobs[0] = (yout[:, species] * factors).sum(1)

    #transpose yobs to make it easy to plot
    yobs.T
    return (xout,yobs,yout)

def read_csv_array(fp):
    """returns the first string and a numpy array from a csv set of data"""
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
    #assume all entries in the list are the same length
    darray = numpy.zeros((len(templist[0]), len(templist)))
    for i, item in enumerate(templist):
        darray[i] = numpy.asarray(templist[i], dtype=darray.dtype)    
    return (headstring, darray)


def compare_data(array0, array1, array0var=None, obspec=None):
    """Compares two arrays of different size and returns the X^2 between them.
    Uses the X axis as the unit to re-grid both arrays
    obspec: passed by odeanneal to know which observable to use for comparison
    """
    # figure out the array shapes
    # this expects arrays of the form array([time, measurements])
    # the time is assumed to be roughly the same for both and the 
    # shortest time will be taken as reference to regrid the data
    # the regridding is done using a b-spline interpolation
    # array0var shuold be the variances at every time point
    #

    # sanity checks
    # make sure we are comparing the right shape arrays
    arr0shape = array0.shape
    arr1shape = array1.shape
    
    if len(arr0shape) != len(arr1shape):
        raise SystemExit("comparing arrays of different dimensions")
    
    # get the time range where the arrays overlap
    rngmin = max(array0[:,0].min(), array1[:,0].min)
    rngmax = min(array0[:,0].max(), array1[:,0].max)
    rngmin = round(rngmin, -1)
    rngmax = round(rngmax, -1)
    
    # use the experimental gridpoints from the reference array as
    # the new gridset for the model array. notice the time range
    # of the experiment has to be within the model range
    #
    iparray = numpy.zeros(array0.shape)
    iparray[:,0] =  array0[:,0]
    
    # now create a b-spline of the data and fit it to desired range
    tck = scipy.interpolate.splrep(array1[:,0], array1[:,1])
    iparray[:,1] = scipy.interpolate.splev(iparray[:,0], tck)
    
    # we now have x and y axis for the points in the model array
    # calculate the objective function
    
    diffarray = array0[:,1] - iparray[:,1]
    diffsqarray = diffarray * diffarray
    
    # assume a default .05 variance
    if array0var is None:
        array0var = numpy.ones(iparray[:,1].shape)
        array0var = array0var*.05
    
    array0var = numpy.multiply(array0var,array0var)
    array0var = array0var*0.5

    objarray = diffsqarray * array0var
    
    return objarray.sum()
    
