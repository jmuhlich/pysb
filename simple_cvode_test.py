# Skeleton to use pysundials CVODES
#
# 1. Import the necessary stuff
from pysundials import cvodes #CVODES
import nvecserial             #NVECSERIAL for easy NVector manimulation
import ctypes                 #for convenience

# 2. Define structure to hold parameter sensitivities and optional data
class UserData(ctypes.Structure):
    _fields_ = [
        ('p', cvodes.realtype*4)
        ]

PUserData = ctypes.POINTER(UserData)

# 3. Define right-hand side function
# must contain four parameters:
# (t, Nvector_curr, Nvector_new, 

