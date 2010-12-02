from pysundials import cvode
import ctypes
from matplotlib import pyplot

V1 = 1.0
V2 = 10.0
V3 = 1.0
V4 = 1.0

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

S1 = 0
S2 = 1
#S3 = 2

def R1(y):
    return (V1 /( K1_X0 * K1_S2 )*( X0*y[S2 ]))/ \
        (1 + X0/ K1_X0 + y[S2 ]/ K1_S2 +\
             (X0*y[S2 ])/( K1_X0 * K1_S2 ))

def R2(y):
    return (V2 /( K2_S3 * K2_X6 ))*(( icsum -y[S2 ])* X6 )/\
        (1 + (icsum -y[S2 ])/ K2_S3 + X6/ K2_X6 \
             + ((( icsum -y[S2 ])* X6 )/( K2_S3 * K2_X6 )))
 
def R3(y):
    return (V3*y[S1 ])/( y[S1] + K3_S1 )

def R4(y):
    return (V4*y[S1 ])/( y[S1] + K4_S1 )
 
def f(t, y, ydot , f_data ):
    ydot [S2] = R2(y) - R1(y)
    ydot [S1] = R1(y) - R3(y) - R4(y)
    return 0

y = cvode.NVector ([1.0 , 0.7])
icsum = 0.7 + 0.3 # setting the initial sum of the conservation

cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SS, 1.0e-8, 1.0e-12)
cvode.CVDense (cvode_mem , 2)

results = ([] , [], [], [])

t = cvode.realtype(0.0)
results [0]. append (t.value)
results [1]. append (y[S1])
results [2]. append (y[S2])
results [3]. append ((icsum- y[S2]))

iout = 1
tout = 0.1
while iout <= 50:
    ret = cvode.CVode(
        cvode_mem,
        tout,
        y,
        ctypes.byref(t),
        cvode.CV_NORMAL
        )
    if ret != 0:
        print " CVODE Error: %i"%(ret)
        break
    results[0].append(tout)
    results[1].append(y[S1])
    results[2].append(y[S2])
    results[3].append((icsum -y[S2]))
     
    iout += 1
    tout += 0.1

pyplot.title("PySUNDIALS/CVODE")
pyplot.xlabel("t")
pyplot.ylabel("Concentration")
pyplot.plot(
    results[0], results[1], "k-",
    results[0], results[2], "k--",
    results[0], results[3], "k-."
    )
pyplot.legend(("$S_1$", "$S_2$", "$S_3$"))
pyplot.show()
