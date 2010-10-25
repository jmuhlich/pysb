from pysundials import cvode
import ctypes
import nvecserial
from matplotlib import pyplot

# Simple reaction:
#       k1f    k1c
# E + S <-> ES -> E + P
#       k1r

# Parameters
k1f = 1.0e-6  #units: of mol^-1 s^-1
k1r = 1.0e-3  #units: s^-1
k2  = 1.0e-2  #units: s^-1

# Initial specs
E0  = 10000 # y[0] init
S0  = 20000 # y[1] init
ES0 = 0     # y[2] init
P0  = 0     # y[3] init

def Edot(y):
    return (-k1f * y[0] * y[1]) + (k1r * y[2]) + (k2 * y[2])

def Sdot(y):
    return (-k1f * y[0] * y[1]) + (k1r * y[2])

def ESdot(y):
    return (k1f * y[0] * y[1]) - (k1r * y[2]) - (k2 *y[2])

def Pdot(y):
    return (k2 * y[2])

# Right hand side definition:
# The function specified takes the following parameters
# time_step (float)       the current value of the independent variable
# y (NVector)             the vector of current dependent values
# ydot (NVector)          undefined values, contents should be set to the new values of y
# f_data (c_void_p)       pointer to user data set by CVodeSetFdata
def f(t, y, ydot, f_data):
    ydot[0] = Edot(y)
    ydot[1] = Sdot(y)
    ydot[2] = ESdot(y)
    ydot[3] = Pdot(y)
    return 0

y = cvode.NVector([E0, S0, ES0, P0]) #initial values for integration

#necessary steps for SUNDIALS
# Sets up a CVODE internal integrator structure, and returns a handle it in the form of a CVodeMemObj
# CVodeCreate takes two arguments, lmm (CV_BDF or CV_ADAMS) and iter(CV_NEWTON or CV_FUNCTIONAL)
# lmm:   The user of the CVODE package specifies whether to use the
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
#        CVODE linear solver. CV_NEWTON is recommended in case of
#        stiff problems.
cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON) #stiff setup

# CVodeMalloc allocates and initializes memory for a problem to be solved by CVODE
# CVodeMalloc(cvodememobj, func, t0, y0, itol, reltol, abstol)
# cvodememobj: as returned by CVodeCreate()
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
cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, 1.0e-8, 1.0e-12)

#
# call to the CVDense function links the main CVODE integrator with the CVDENSE linear solver
# CVDense(cvodememobj, N)
# cvodememobj: the memobject from CVodeCreate()
# N: (int) the size (number of DE's?) of the ODE system
cvode.CVDense(cvode_mem, 4)

output = ([], [], [], [], [])

t = cvode.realtype(0.0)
results = ([], [], [], [], [])
results[0].append(t.value)
results[1].append(y[0])
results[2].append(y[1])
results[3].append(y[2])
results[4].append(y[3])

iout = 0
tout = 0.1
print "Beginning integration..."
while iout < 10000:
# call to the CVODE integrator
# CVode(cvodememobj, tout, yout, tret, itask)
# cvodememobj: from CVodeCreate()
# tout: (float), next time at which solution is desired
# yout: (NVector) computed solution vector. in CV_NORMAL with no errors or roots found
# tret: (*realtype) is set to the time reached by the solver
# itask: (int) is one of CV_NORMAL, CV_ONE_STEP, CV_NORMAL_TSTOP, or CV_ONE_STEP_TSTOP
#
    ret = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)

    if ret != 0:
        print "CVODE ERROR: %i"%(ret)
        break

    results[0].append(tout)
    results[1].append(y[0])
    results[2].append(y[1])
    results[3].append(y[2])
    results[4].append(y[3])
    
    iout += 1
    tout += 0.1
print "Integration finished."

pyplot.plot(
    results[0], results[1], 'k-',
    results[0], results[2], 'k--',
    results[0], results[3], 'k-.',
    results[0], results[4], 'k:'
    )
pyplot.legend(('E', 'S', 'ES', 'P'))

    
