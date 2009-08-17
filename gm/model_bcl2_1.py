from pysb import *


Model('model')


Parameter('k_baxm_dim_f', 1)
Parameter('k_baxm_dim_r', 1)
Parameter('k_bak_dim_f',  1)
Parameter('k_bak_dim_r',  1)


Monomer('CASP8',  ['sub', 'mcl1', 'active'],         {'active': ['no', 'yes']})
Monomer('BIDC',   ['csp8', 'pro', 'anti', 'active'], {'active': ['no', 'yes']})
Monomer('BIDM',   ['csp8', 'pro', 'anti', 'active'], {'active': ['no', 'yes']})
Monomer('BAXC',   ['activator', 'anti', 'active'],   {'active': ['no', 'yes']})
Monomer('BAXM',   ['dim', 'tetra', 'activator', 'anti', 'active'],
                                                     {'active': ['no', 'yes']})
Monomer('BIM',    ['pro', 'anti'])
Monomer('BAK',    ['dim', 'tetra', 'activator', 'anti', 'active'],
                                                     {'active': ['no', 'yes']})
Monomer('BCLXLC', ['bh3', 'pro'])
Monomer('BCLXLM', ['bh3', 'pro'])
Monomer('BCL2',   ['bh3', 'pro'])


Rule('BAXM_dimerize',
     BAXM(dim=None, active='yes') + BAXM(dim=None, active='yes') <> BAXM(dim=1, active='yes') * BAXM(dim=1, active='yes'),
     k_baxm_dim_f, k_baxm_dim_r)

Rule('BAK_dimerize',
     BAK(dim=None, active='yes') + BAK(dim=None, active='yes') <> BAK(dim=1, active='yes') * BAK(dim=1, active='yes'),
     k_bak_dim_f, k_bak_dim_r)


model.observe('baxm_dimer', BAXM(active='yes') * BAXM(active='yes'))
model.observe('bak_dimer', BAK(active='yes') * BAK(active='yes'))
