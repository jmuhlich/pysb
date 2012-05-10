from pysb import *

# from http://bionetgen.org/index.php/Simple

Model()


# Physical and geometric constants
Parameter('NA', 6.0e23)     # Avogadro's num
Parameter('f', 0.01)        # scaling factor
Parameter('Vo', f * 1e-10)  # L
Parameter('V', f * 3e-12)   # L
 
# Initial concentrations
Parameter('EGF_conc', 2e-9)            # nM
Parameter('EGF0', EGF_conc * NA * Vo)  # nM
Parameter('EGFR0', f * 1.8e5)          # copy per cell

# Rate constants 
Parameter('kp1', 9.0e7 / (NA * Vo))  # input /M/sec         
Parameter('km1', 0.06)               # /sec    


Monomer('EGF', ['R'])
Monomer('EGFR', ['L', 'CR1', 'Y1068'], {'Y1068': ['U', 'P']})


Initial(EGF(R=None), EGF0)
Initial(EGFR(L=None, CR1=None, Y1068='U'), EGFR0)


Rule('egf_binds_egfr', EGF(R=None) + EGFR(L=None) <> EGF(R=1) % EGFR(L=1), kp1, km1)


# Species LR EGF(R!1).EGFR(L!1)
Observable('Lbound', EGF(R=ANY))  # Molecules


if __name__ == '__main__':
    from pylab import linspace, plot, legend, show
    from pysb.integrate import odesolve
    t = linspace(0, 100, 20)
    yrec = odesolve(model, t)
    y = yrec.view().reshape(len(yrec), len(yrec.dtype))
    plot(t, y)
    legend(model.species + model.observables.values())
    show()
