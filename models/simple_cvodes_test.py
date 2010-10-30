from pysundials import cvodes, nvecserial
import ctypes
from matplotlib import pyplot

# Simple reaction:
#       k1f    k1c
# E + S <-> ES -> E + P
#       k1r

# structure to hold sensitivity data parameters
class UserData(ctypes.Structure):
    _fields_ = [('p', cvodes.realtype*3)] #3 parameters

PUserData = ctypes.POINTER(UserData)

# Parameters
#k1f = 1.0e-6  #units: of mol^-1 s^-1
#k1r = 1.0e-3  #units: s^-1
#k2  = 1.0e-2  #units: s^-1
k1f = 0
k1r = 1
k2  = 2

data = UserData() 
data.p[k1f] = 1.0e-6
data.p[k1r] = 1.0e-3
data.p[k2]  = 1.0e-2

# Initial specs
E0  = 10000 # y[0] init
S0  = 20000 # y[1] init
ES0 = 0     # y[2] init
P0  = 0     # y[3] init

def Edot(y, p):
    return (-p[k1f] * y[0] * y[1]) + (p[k1r] * y[2]) + (p[k2] * y[2])

def Sdot(y, p):
    return (-p[k1f] * y[0] * y[1]) + (p[k1r] * y[2])

def ESdot(y, p):
    return (p[k1f] * y[0] * y[1]) - (p[k1r] * y[2]) - (p[k2] *y[2])

def Pdot(y, p):
    return (p[k2] * y[2])

edot = 0
sdot = 1
esdot = 2
pdot = 3

# Right hand side definition:
# The function specified takes the following parameters
# time_step (float)       the current value of the independent variable
# y (NVector)             the vector of current dependent values
# ydot (NVector)          undefined values, contents should be set to the new values of y
# f_data (c_void_p)       pointer to user data set by CvodesSetFdata
def f(t, y, ydot, f_data):
    data = ctypes.cast(f_data, PUserData).contents
    ydot[edot] = Edot(y, data.p)
    ydot[sdot] = Sdot(y, data.p)
    ydot[esdot] = ESdot(y, data.p)
    ydot[pdot] = Pdot(y, data.p)
    return 0

y = cvodes.NVector([E0, S0, ES0, P0]) #initial values for integration

#necessary steps for SUNDIALS
# Sets up a CVODES internal integrator structure, and returns a handle it in the form of a CvodesMemObj
# CvodesCreate takes two arguments, lmm (CV_BDF or CV_ADAMS) and iter(CV_NEWTON or CV_FUNCTIONAL)
# lmm:   The user of the CVODES package specifies whether to use the
#        CV_ADAMS (Adams-Moulton) or CV_BDF (Backward Differentiation
#        Formula) linear multistep method. The BDF method is
#        recommended for stiff problems, and the CV_ADAMS method is
#        recommended for nonstiff problems.
#
# iter:  At each internal time step, a nonlinear equation must
#        be solved. The user can specify either CV_FUNCTIONAL
#        iteration, which does not require linear algebra, or a
#        CV_NEWTON iteration, which requires the solution of linear
#        systems. In the CV_NEWTON case, the user also specifies a
#        CVODES linear solver. CV_NEWTON is recommended in case of
#        stiff problems.
cvodes_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON) #stiff setup

# CvodesMalloc allocates and initializes memory for a problem to be solved by CVODES
# CvodesMalloc(cvodesmemobj, func, t0, y0, itol, reltol, abstol)
# cvodesmemobj: as returned by CvodesCreate()
# func: python callable function defining right hand side of fxn y' = f(t,y)
# t0: initial value of dependent variable (time in this case)
# y0: initial value of condition vector y(t0)
# itol: type of tolerances method to be used:
#       CV_SS - scalar relative and absolute tolerances
#       CV_SV - scalar relative tolerance and a vector of absolute tolerances.
#       CV_WF - indicates that the user will provide a function to evaluate the error weights. 
#               In this case, reltol and abstol are ignored, and should be None.
# reltol: (float) the relative tolerance scalar
# abstol: (float/NVector) absolute tolerance(s) 
#
#abstol = cvodes.NVector([1.0e-8, 1.0e-8, 1.0e-6])
#abstol = cvodes.NVector([1.0e-2, 1.0e-2, 1.0e-2])
abstol = 1.0e-12
#reltol = cvodes.realtype(1.0e-4)
reltol = 1.0e-8
cvodes.CVodeMalloc(cvodes_mem, f, 0.0, y, cvodes.CV_SS, reltol, abstol)

# Set optional inputs (in this case the sensitivity data pointers)
#
cvodes.CVodeSetFdata(cvodes_mem, ctypes.pointer(data)) #points to sensitivity data

# call to the CVDense function links the main CVODES integrator with the CVDENSE linear solver
# CVDense(cvodesmemobj, N)
# cvodesmemobj: the memobject from CvodesCreate()
# N: (int) the size (number of DE's?) of the ODE system
cvodes.CVDense(cvodes_mem, 4)

# set sensitivity system options 
yS = nvecserial.NVectorArray([([0]*4)]*3)

print 'yS:', yS

# CVodeSensMalloc(cvodememobj, Ns, ism, yS0)
#    CVodeSensMalloc allocates and initializes memory related to sensitivity computations.
#    cvodememobj: cvodes mem object 
#    Ns (int): number of sensitivities to be computed
#    ism (int):type of corrector used in sensitivity analysis. 
#    yS0 (NVectorArray)              is the array of initial condition vectors for sensitivity variables

cvodes.CVodeSensMalloc(cvodes_mem, 3, cvodes.CV_SIMULTANEOUS, yS)

# CVodeSetSensParams expects four parameters
# (for more detail see p. 111 of the CVODES user guide)
# 1. the cvodes memory object
# 2. a pointer to the array of parameter values which MUST be passed
#    through the user data structure (so CVODES knows where the values
#    are and can peturb them, presumably)
# 3. an array (i.e. list) of scaling factors, one for each parameter
#    for which sensitivies are to be determined
# 4. an array of integers (either 1 or 0) , where a 1 indicates the
#    respective parameter value should be used in estimating
#    sensitivities

cvodes.CVodeSetSensParams(cvodes_mem,
                          data.p, 
                          [1, 1, 1],
                          [1, 1, 1]
                          )
t = cvodes.realtype(0.0)
results = [[] for i in range(0,17)]
results[0].append(t.value)
results[1].append(y[edot])
results[2].append(y[sdot])
results[3].append(y[esdot])
results[4].append(y[pdot])
results[5].append(yS[k1f][edot])  #dE/dk1f 
results[6].append(yS[k1r][edot])  #dE/dk1r 
results[7].append(yS[k2][edot])   #dE/dk2 
results[8].append(yS[k1f][sdot])  #dS/dk1f 
results[9].append(yS[k1r][sdot])  #dS/dk1r 
results[10].append(yS[k2][sdot])  #dS/dk2 
results[11].append(yS[k1f][esdot]) #dES/dk1f 
results[12].append(yS[k1r][esdot]) #dES/dk1r 
results[13].append(yS[k2][esdot]) #dES/dk2 
results[14].append(yS[k1f][pdot])  #dP/dk1f 
results[15].append(yS[k1r][pdot])  #dP/dk1r 
results[16].append(yS[k2][pdot])  #dP/dk2 

iout = 0
tout = 0.1
print "Beginning integration..."
while iout < 10000:
# call to the CVODES integrator
# Cvodess(cvodesmemobj, tout, yout, tret, itask)
# cvodesmemobj: from CvodesCreate()
# tout: (float), next time at which solution is desired
# yout: (NVector) computed solution vector. in CV_NORMAL with no errors or roots found
# tret: (*realtype) is set to the time reached by the solver
# itask: (int) is one of CV_NORMAL, CV_ONE_STEP, CV_NORMAL_TSTOP, or CV_ONE_STEP_TSTOP
#
    #print 'yS:\n', yS

    ret = cvodes.CVode(cvodes_mem, tout, y, ctypes.byref(t), cvodes.CV_NORMAL)
    cvodes.CVodeGetSens(cvodes_mem, t, yS)
    

    if ret != 0:
        print "CVODES ERROR: %i"%(ret)
        break

    results[0].append(tout)
    results[1].append(y[edot])
    results[2].append(y[sdot])
    results[3].append(y[esdot])
    results[4].append(y[pdot])
    results[5].append(yS[k1f][edot])  #dE/dk1f 
    results[6].append(yS[k1r][edot])  #dE/dk1r 
    results[7].append(yS[k2][edot])   #dE/dk2 
    results[8].append(yS[k1f][sdot])  #dS/dk1f 
    results[9].append(yS[k1r][sdot])  #dS/dk1r 
    results[10].append(yS[k2][sdot])  #dS/dk2 
    results[11].append(yS[k1f][esdot])#dES/dk1f 
    results[12].append(yS[k1r][esdot])#dES/dk1r 
    results[13].append(yS[k2][esdot]) #dES/dk2 
    results[14].append(yS[k1f][pdot]) #dP/dk1f 
    results[15].append(yS[k1r][pdot]) #dP/dk1r 
    results[16].append(yS[k2][pdot])  #dP/dk2 
    
    iout += 1
    tout += 0.1
print "Integration finished."

pyplot.figure(1)
pyplot.subplot(511)
pyplot.ylabel("Concentration")
pyplot.plot(
    results[0], results[1], 'k-',
    results[0], results[2], 'k--',
    results[0], results[3], 'k-.',
    results[0], results[4], 'k:'   
    )
pyplot.legend(('E', 'S', 'ES', 'P'))

pyplot.subplot(512)
pyplot.plot(
    results[0], results[5], 'k-',  
    results[0], results[6], 'k--', 
    results[0], results[7], 'k-.'
    )
pyplot.legend(('$\delta E/dk1f $', '$\delta E/dk1r $', '$\delta E/dk2 $'))

pyplot.subplot(513)
pyplot.plot(
    results[0], results[8], 'k-',  
    results[0], results[9], 'k--', 
    results[0], results[10], 'k-.'
    )
pyplot.legend(('$\delta S/dk1f $', '$\delta S/dk1r $', '$\delta S/dk2 $'))

pyplot.subplot(514)
pyplot.plot(
    results[0], results[11], 'k-',  
    results[0], results[12], 'k--', 
    results[0], results[13], 'k-.'
    )
pyplot.legend(('$\delta ES/dk1f $', '$\delta ES/dk1r $', '$\delta ES/dk2 $'))

pyplot.subplot(515)
pyplot.xlabel("t")
pyplot.plot(
    results[0], results[14], 'k-',  
    results[0], results[15], 'k--', 
    results[0], results[16], 'k-.'
    )
pyplot.legend(('$\delta P/dk1f $', '$\delta P/dk1r $', '$\delta P/dk2 $'))

#pyplot.show()
