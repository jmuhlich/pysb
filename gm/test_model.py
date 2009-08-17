import sys

name = 'model_bcl2_1'
if len(sys.argv) > 1:
    name = sys.argv[1]
m = __import__(name)

import pysb.generator.bng as bng
gen = bng.BngGenerator(m.model)
print gen.content
