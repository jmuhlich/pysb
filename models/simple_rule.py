from pysb import *

#simple test for E + S <-> ES -> E + P

Model()

Parameter('k1f', 5.0e-5)
Parameter('k1r', 1.0e-3)
Parameter('k1c', 1.0e-1)

Monomer("Enz", 'b')
Monomer("Sub", ['b','state'], {'state':['i','a']})

Rule('ES_binding',
     Enz(b=None) + Sub(b=None, state='i') <>
     Enz(b=1   ) % Sub(b=1   , state='i'),
     k1f, k1r)

Rule('Sub_activation',
     Enz(b=1   ) % Sub(b=1   , state='i') >>
     Enz(b=None) + Sub(b=None, state='a'),
     k1c)

#initial amounts
Parameter('Enz_0', 10000); # Initial enzyme
Parameter('Sub_0', 10000); # Initial substrate
Initial(Enz(b=None), Enz_0)
Initial(Sub(b=None, state='i'), Sub_0)

#Observables
Observe('Enz', Enz(b=None))
Observe('Sub', Sub(b=None, state='i'))
Observe('EnzSub', Enz(b=1)%Sub(b=1))
Observe('Prod', Sub(b=None, state = 'a'))

# # generate initial conditions from _0 parameter naming convention
# # for m in model.monomers:
# #     ic_param = model.parameter('%s_0' % m.name)
# #     if ic_param is not None:
# #         sites = {}
# #         for s in m.sites:
# #             if s in m.site_states:
# #                 sites[s] = m.site_states[s][0]
# #             else:
# #                 sites[s] = None
# #         Initial(m(sites), ic_param)

# if __name__ == '__main__':
#     from pysb.bng import generate_network_code
#     from pysb.tools.export_bng import run as run_export
#     print run_export(model)

