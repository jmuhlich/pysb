from pysundials import cvode
import ctypes

# create an Nvector

a = cvode.NVector([1, 2, 3])

print "this is the NVector:", a

print "this is the NVector + 1", a+1


