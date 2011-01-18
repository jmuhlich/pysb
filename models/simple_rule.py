from pysb import *

#simple test for E + S <-> ES -> E + P

Model()
#initial amounts
Parameter('Enz_0', 10000); # Initial enzyme
Parameter('Sub_0', 20000); # Initial substrate

Parameter('k1f', 1.0e-6)
Parameter('k1r', 1.0e-3)
Parameter('k1c', 1.0e-2)

Monomer("Enz", 'b')
Monomer("Sub", ['b','state'], {'state':['i','a']})

Rule('ES binding',
     Enz(b=None) + Sub(b=None, state='i') <>
     Enz(b=1   ) % Sub(b=1   , state='i'),
     k1f, k1r)

Rule('Sub activation',
     Enz(b=1   ) % Sub(b=1   , state='i') >>
     Enz(b=None) % Sub(b=None, state='a'),
     k1c)


Observe('Enz', Enz(b=None))
Observe('Sub', Sub(b=None))
Observe('EnzSub', Enz(b=1)%Sub(b=1, state='i'))
Observe('Prod', Sub(b=None))

# generate initial conditions from _0 parameter naming convention
for m in model.monomers:
    ic_param = model.parameter('%s_0' % m.name)
    if ic_param is not None:
        sites = {}
        for s in m.sites:
            if s in m.site_states:
                sites[s] = m.site_states[s][0]
            else:
                sites[s] = None
        Initial(m(sites), ic_param)
Initial(Enz(b=None), Enz_0)
Initial(Sub(b=None, state='i'), Sub_0)

if __name__ == '__main__':
    from pysb.generator.bng import BngGenerator
    gen = BngGenerator(model)
    print gen.get_content()
    print ""
    print "begin actions"
    print "  generate_network({overwrite=>1});"
    print "  simulate_ode({t_end=>21600,n_steps=>360});" # 6 hours, 1-minute steps
    print "end actions"
