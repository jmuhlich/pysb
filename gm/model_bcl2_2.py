from pysb import *


def can_activate(*monomers):
    for m in monomers:
        m.sites += ['active']
        m.sites_dict['active'] = None
        m.site_states['active'] = ['no', 'yes']

def dimerize(monomer, rate_f, rate_r):
    return Rule(monomer.name + '_dimerize',
                monomer(dim=None, active='yes') + monomer(dim=None, active='yes') <>
                monomer(dim=1,    active='yes') * monomer(dim=1,    active='yes'),
                rate_f, rate_r)


Model('model')


Parameter('k_baxm_dim_f', 1)
Parameter('k_baxm_dim_r', 1)
Parameter('k_bak_dim_f',  1)
Parameter('k_bak_dim_r',  1)


Monomer('CASP8',  ['sub', 'mcl1'])
Monomer('BIDC',   ['csp8', 'pro', 'anti'])
Monomer('BIDM',   ['csp8', 'pro', 'anti'])
Monomer('BAXC',   ['activator', 'anti'])
Monomer('BAXM',   ['dim', 'tetra', 'activator', 'anti'])
Monomer('BIM',    ['pro', 'anti'])
Monomer('BAK',    ['dim', 'tetra', 'activator', 'anti'])
Monomer('BCLXLC', ['bh3', 'pro'])
Monomer('BCLXLM', ['bh3', 'pro'])
Monomer('BCL2',   ['bh3', 'pro'])

can_activate(CASP8, BIDC, BIDM, BAXC, BAXM, BAK)


dimerize(BAXM, k_baxm_dim_f, k_baxm_dim_r)
dimerize(BAK,  k_bak_dim_f,  k_bak_dim_r)


model.observe('baxm_dimer', BAXM(active='yes') * BAXM(active='yes'))
model.observe('bak_dimer', BAK(active='yes') * BAK(active='yes'))
