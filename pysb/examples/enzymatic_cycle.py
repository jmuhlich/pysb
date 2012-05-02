from pysb import *

# from http://bionetgen.org/index.php/Enzymatic_Cycle

Model()


# kinetic parameters
Parameter('k1', 1.0)    # units: /time/molecule
Parameter('k2', 1.0)    # units: /time
Parameter('k3', 0.1)    # units: /time
Parameter('k4', 1.0)    # units: /time/molecule
Parameter('k5', 1.0)    # units: /time
Parameter('k6', 0.1)    # units: /time
# total enzyme and substrate (units: molecules)
Parameter('E1_tot', 1)
Parameter('E4_tot', 1)
Parameter('S2_tot', 50)
Parameter('S5_tot', 50)
# derived parameters (units: molecules)
Parameter('Km2', (k2 + k3)/k1)
Parameter('Km5', (k5 + k6)/k4)
# calculate initial free S2 and S5, assuming quasi-equillibrium enzyme binding
Parameter('S2_free0', (S2_tot + E1_tot - Km2)/2 + ((S2_tot + E1_tot - Km2) ** 2 + 4*Km2)**0.5/2)
Parameter('S5_free0', (S5_tot + E4_tot - Km5)/2 + ((S5_tot + E4_tot - Km5) ** 2 + 4*Km5)**0.5/2)


Monomer('S', ['c'], {'c': ['c2', 'c5']})
Monomer('E1')
Monomer('E4')


Initial(E1(), E1_tot)           # S1+S3 (total enzyme A)
Initial(S(c='c2'), S2_free0)    # S2 (free substrate in config 2)
Initial(E4(), E4_tot)           # S4+S6 (total enzyme B)
Initial(S(c='c5'), S5_free0)    # S5 (free substrate in config 5)


# free substrate (required for functional ratelaw)
Observe('S2_free', S(c='c2'))
Observe('S5_free', S(c='c5'))


# TODO: make these functional rate laws work
# s2 -> s5, catalyzed by E1
Rule('s2_s5_cat_E1', S(c='c2') + E1() >> S(c='c5') + E1(), k3/(Km2 + S2_free))
# s5 -> s2, catalyzed by E4
Rule('s5_s2_cat_E4', S(c='c5') + E4() >> S(c='c2') + E4(), k6/(Km5 + S5_free))
