from pysundials import cvode
import ctypes
from matplotlib import pyplot

# Simple reaction:
# E + S <-> ES -> E + P

# Parameters
k1f = 2    #units: of mol^-1 s^-1
k1r = 1    #units: s^-1
k1r = 1.5  #units: s^-1

# Initial specs
S  = 80000
E  = 40000
ES = 0
P  = 0

def R1F(y):
    return 
