from pysundials import cvodes
from pysundials import nvecserial
import ctypes
from matplotlib import pyplot

K1_X0 = 1.0
K1_S2 = 1.0
K2_S3 = 1.0
K2_X6 = 1.0
K3_S1 = 1.0
K4_S1 = 1.0

X0 = 1.0
X4 = 1.0e-12
X5 = 1.0e-12
X6 = 1.0
X7 = 1.0e-12

#indexes for substances
S1 = 0
S2 = 1

#these are now indexes, not real values
V1 = 0
V2 = 1
V3 = 2
V4 = 3

class UserData (ctypes.Structure):
    _fields_ = [("p", cvodes.realtype*4)]

PUserData = ctypes.POINTER(UserData)

def R1(y, p):
    return (p[V1] /( K1_X0 * K1_S2 )*( X0*y[S2])) / \
        (1 + X0/ K1_X0 + y[S2]/ K1_S2 + (X0*y[S2])/( K1_X0 * K1_S2 ))

def R2(y, p):
    return (p[V2] /( K2_S3 * K2_X6 ))*(( icsum -y[S2])* X6 )/\
        (1 + (icsum -y[S2])/ K2_S3 + X6/ K2_X6 \
             + ((( icsum - y[S2])* X6)/( K2_S3 * K2_X6 )))
 
def R3(y, p):
    return (p[V3]*y[S1])/(y[S1] + K3_S1)

def R4(y, p):
    return (p[V4]*y[S1])/(y[S1] + K4_S1)
 
def f(t, y, ydot, f_data):
    data = ctypes.cast(f_data, PUserData).contents

    ydot [S2] = R2(y, data.p) - R1(y, data.p)
    ydot [S1] = R1(y, data.p) - R3(y, data.p) - R4(y, data.p)
    return 0
# The previous assignment was:
# V1 = 1.0
# V2 = 10.0
# V3 = 1.0
# V4 = 1.0
# Now they are added into the new data.p object
data = UserData() #define an instance of UserData
data.p[V1] = 1.0
data.p[V2] = 10.0
data.p[V3] = 1.0
data.p[V4] = 1.0

y = cvodes.NVector([1.0 , 0.7])
icsum = 0.7 + 0.3 # setting the initial sum of the conservation

cvodes_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
#
# Set CVodeMalloc(cvode_mem, f, 0.0, y, cvodes.CV_SV, reltol, abstol)
# reltol and abstol can be values or vectors
# if using a vector then cvodes.CV_SV shuold be used
#
cvodes.CVodeMalloc(cvodes_mem, f, 0.0, y, cvodes.CV_SS, 1.0e-8, 1.0e-12)
cvodes.CVodeSetFdata(cvodes_mem, ctypes.pointer(data)) #point to the sens data
cvodes.CVDense(cvodes_mem, 2) #choice of linear solver, cvdense, cvband, etc
 
#
#create NvectorArray of size "# of variables" x "# of params for sensitivities"
#
yS = nvecserial.NVectorArray([([0]*2)]*4)

# 
#
cvodes.CVodeSensMalloc(cvodes_mem, 4, cvodes.CV_SIMULTANEOUS, yS)

# So this bit below might seem like magic , but it's fairly straight
# forward. This particular line doesn't look like the example
# cvsfwddenx.py because we are not suplying our own sensitivity
# function.

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
                          data.p, # we have four system parameters ( The four V's)
                          [1, 1, 1, 1],  # they should all be scaled by 1, i.e. unscaled,
                          [1, 1, 1, 1]   # contrib to sensitivities estimation
                          )

results = ([], [], [], [], [], [], [], [], [], [])

t = cvodes.realtype(0.0)
results[0].append(t.value)
results[1].append(y[S1])
results[2].append(y[S2])
results[3].append((icsum- y[S2]))
results[4].append(yS[V2][S1])
results[5].append(yS[V2][S2])
results[6].append(yS[V3][S1])
results[7].append(yS[V3][S1])
results[8].append(yS[V4][S1])
results[9].append(yS[V4][S2])
iout = 1
tout = 0.05
while iout <= 80:
    ret = cvodes.CVode(cvodes_mem,
                       tout,
                       y,
                       ctypes.byref(t),
                       cvodes.CV_NORMAL
                       )
    cvodes.CVodeGetSens(cvodes_mem, t, yS)
    if ret != 0:
        print " CVODE Error: %i"%(ret)
        break
    results[0].append(tout)
    results[1].append(y[S1])
    results[2].append(y[S2])
    results[3].append((icsum -y[S2]))
    results[4].append(yS[V2][S1])
    results[5].append(yS[V2][S2])
    results[6].append(yS[V3][S1])
    results[7].append(yS[V3][S1])
    results[8].append(yS[V4][S1])
    results[9].append(yS[V4][S2])
    iout += 1
    tout += 0.1

pyplot.figure(1)
pyplot.subplot(411)
pyplot.ylabel("Concentration")
pyplot.plot(
    results[0], results[1], "k-",
    results[0], results[2], "k--",
    results[0], results[3], "k-."
    )
pyplot.legend(("$S_1$", "$S_2$", "$S_3$"))

pyplot.subplot(412)
pyplot.xlabel ("t")
pyplot.plot(
    results[0], results[4], "k-",
    results[0], results[5], "k--",
    )
pyplot.legend(("$\delta S_1/dV_2$", "$\delta S_2/dV_2$"))

pyplot.subplot(413)
pyplot.xlabel ("t")
pyplot.plot(
    results[0], results[4], "k-",
    results[0], results[5], "k--",
    )
pyplot.legend(("$\delta S_1/dV_3$", "$\delta S_2/dV_3$"))

pyplot.subplot(414)
pyplot.xlabel ("t")
pyplot.plot(
    results[0], results[8], "b-",
    results[0], results[9], "b--",
    )
pyplot.legend(("$\delta S_1/dV_4$", "$\delta S_2/dV_4$"))

pyplot.show()
