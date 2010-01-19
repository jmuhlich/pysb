from pysb import *

# Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK, 2008
# Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic
# Cell Death. PLoS Biol 6(12): e299. doi:10.1371/journal.pbio.0060299
#
# http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0060299

Model()

transloc = .01; # rate of translocation between the cytosolic and mitochondrial compartments

v = .07; # mitochondria compartment volume/cell volume

# Non-zero initial conditions (in molecules per cell):
Parameter('Lsat'       , 6e4);  # saturating level of ligand (corresponding to ~1000 ng/ml SuperKiller TRAIL)
Parameter('L_0'        , 3000); # baseline level of ligand for most experiments (corresponding to 50 ng/ml SuperKiller TRAIL)
Parameter('pR_0'       , 200);  # TRAIL receptor (for experiments not involving siRNA)
Parameter('pR_0_siRNA' , 100000); # TRAIL receptor for experiments involving siRNA; this is set higher than in non-siRNA experiments to reflect the experimentally observed sensitization caused by both targeting and non-targeting siRNA transfection, which is believed to occur at least partly through the upregulation of death receptors by an interferon pathway.
Parameter('flip_0'     , 1e2);  # Flip
Parameter('pC8_0'      , 2e4);  # procaspase-8 (pro-C8)
Parameter('BAR_0'      , 1e3);  # Bifunctional apoptosis regulator
Parameter('pC3_0'      , 1e4);  # procaspase-3 (pro-C3)
Parameter('pC6_0'      , 1e4);  # procaspase-6 (pro-C6)  
Parameter('XIAP_0'     , 1e5);  # X-linked inhibitor of apoptosis protein  
Parameter('PARP_0'     , 1e6);  # C3* substrate
Parameter('Bid_0'      , 4e4);  # Bid
Parameter('Bcl2c_0'    , 2e4);  # cytosolic Bcl-2
Parameter('Bax_0'      , 1e5);  # Bax
Parameter('Bcl2_0'     , 2e4);  # mitochondrial Bcl-2  
Parameter('Mito_0'     , 5e5);  # mitochondrial binding sites for activated Bax
Parameter('mCytoC_0'   , 5e5);  # cytochrome c
Parameter('mSmac_0'    , 1e5);  # Smac    
Parameter('pC9_0'      , 1e5);  # procaspase-9 (pro-C9)
Parameter('Apaf_0'     , 1e5);  # Apaf-1


Monomer('L', ['b'])
Monomer('pR', ['b'])
Monomer('DISC', ['b'])
Monomer('flip', ['b'])
Monomer('pC8', ['b'])
Monomer('C8', ['b'])
Monomer('BAR', ['b'])
Monomer('pC3', ['b'])
Monomer('C3', ['b'])
Monomer('pC6', ['b'])
Monomer('C6', ['b'])
Monomer('XIAP', ['b'])
Monomer('C3_U', ['b'])
Monomer('PARP', ['b'])
Monomer('CPARP', ['b'])
Monomer('Bid', ['b'])
Monomer('tBid', ['b'])
Monomer('Bcl2c', ['b'])
Monomer('Bax', ['b'])
Monomer('aBax', ['b'])
Monomer('MBax', ['b'])
Monomer('Bcl2', ['b'])
Monomer('Bax2', ['b'])
Monomer('Bax4', ['b'])
Monomer('Mito', ['b'])
Monomer('AMito', ['b'])
Monomer('mCytoC', ['b'])
Monomer('ACytoC', ['b'])
Monomer('mSmac', ['b'])
Monomer('ASmac', ['b'])
Monomer('cCytoC', ['b'])
Monomer('Apaf', ['b'])
Monomer('aApaf', ['b'])
Monomer('pC9', ['b'])
Monomer('Apop', ['b'])
Monomer('cSmac', ['b'])


# L + pR <--> L:pR --> DISC
Parameter('kf1', 4e-07)
Parameter('kr1', 1e-03)
Parameter('kc1', 1e-05)
Rule('bind_L_pR', L(b=None) + pR(b=None) <> L(b=1) * pR(b=1), kf1, kr1)
Rule('produce_DISC', L(b=1) * pR(b=1) >> DISC(b=None), kc1)

# flip + DISC <-->  flip:DISC  
Parameter('kf2', 1e-06)
Parameter('kr2', 1e-03)
Rule('bind_flip_DISC', flip(b=None) + DISC(b=None) <> flip(b=1) * DISC(b=1), kf2, kr2)

# pC8 + DISC <--> DISC:pC8 --> C8 + DISC
Parameter('kf3', 1e-06)
Parameter('kr3', 1e-03)
Parameter('kc3', 1e+00)
Rule('bind_pC8_DISC', pC8(b=None) + DISC(b=None) <> DISC(b=1) * pC8(b=1), kf3, kr3)
Rule('cat_C8_DISC', DISC(b=1) * pC8(b=1) >> C8(b=None) + DISC(b=None), kc3)

# C8 + BAR <--> BAR:C8 
Parameter('kf4', 1e-06)
Parameter('kr4', 1e-03)
Rule('bind_C8_BAR', C8(b=None) + BAR(b=None) <> BAR(b=1) * C8(b=1), kf4, kr4)

# pC3 + C8 <--> pC3:C8 --> C3 + C8
Parameter('kf5', 1e-07)
Parameter('kr5', 1e-03)
Parameter('kc5', 1e+00)
Rule('bind_pC3_C8', pC3(b=None) + C8(b=None) <> pC3(b=1) * C8(b=1), kf5, kr5)
Rule('cat_C3_C8', pC3(b=1) * C8(b=1) >> C3(b=None) + C8(b=None), kc5)

# pC6 + C3 <--> pC6:C3 --> C6 + C3
Parameter('kf6', 1e-06)
Parameter('kr6', 1e-03)
Parameter('kc6', 1e+00)
Rule('bind_pC6_C3', pC6(b=None) + C3(b=None) <> pC6(b=1) * C3(b=1), kf6, kr6)
Rule('cat_C6_C3', pC6(b=1) * C3(b=1) >> C6(b=None) + C3(b=None), kc6)

# pC8 + C6 <--> pC8:C6 --> C8 + C6
Parameter('kf7', 3e-08)
Parameter('kr7', 1e-03)
Parameter('kc7', 1e+00)
Rule('bind_pC8_C6', pC8(b=None) + C6(b=None) <> pC8(b=1) * C6(b=1), kf7, kr7)
Rule('cat_C8_C6', pC8(b=1) * C6(b=1) >> C8(b=None) + C6(b=None), kc7)

# XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U
Parameter('kf8', 2e-06)
Parameter('kr8', 1e-03)
Parameter('kc8', 1e-01)
Rule('bind_XIAP_C3', XIAP(b=None) + C3(b=None) <> XIAP(b=1) * C3(b=1), kf8, kr8)
Rule('cat_XIAP_C3_U', XIAP(b=1) * C3(b=1) >> XIAP(b=None) + C3_U(b=None), kc8)

# PARP + C3 <--> PARP:C3 --> CPARP + C3
Parameter('kf9', 1e-06)
Parameter('kr9', 1e-02)
Parameter('kc9', 1e+00)
Rule('bind_PARP_C3', PARP(b=None) + C3(b=None) <> PARP(b=1) * C3(b=1), kf9, kr9)
Rule('cat_CPARP_C3', PARP(b=1) * C3(b=1) >> CPARP(b=None) + C3(b=None), kc9)

# Bid + C8 <--> Bid:C8 --> tBid + C8
Parameter('kf10', 1e-07)
Parameter('kr10', 1e-03)
Parameter('kc10', 1e+00)
Rule('bind_Bid_C8', Bid(b=None) + C8(b=None) <> Bid(b=1) * C8(b=1), kf10, kr10)
Rule('cat_tBid_C8', Bid(b=1) * C8(b=1) >> tBid(b=None) + C8(b=None), kc10)

# tBid + Bcl2c <-->  tBid:Bcl2c  
Parameter('kf11', 1e-06)
Parameter('kr11', 1e-03)
Rule('bind_tBid_Bcl2c', tBid(b=None) + Bcl2c(b=None) <> tBid(b=1) * Bcl2c(b=1), kf11, kr11)

# Bax + tBid <--> Bax:tBid --> aBax + tBid 
Parameter('kf12', 1e-07)
Parameter('kr12', 1e-03)
Parameter('kc12', 1e+00)
Rule('bind_Bax_tBid', Bax(b=None) + tBid(b=None) <> Bax(b=1) * tBid(b=1), kf12, kr12)
Rule('cat_aBax_tBid', Bax(b=1) * tBid(b=1) >> aBax(b=None) + tBid(b=None), kc12)

# aBax <-->  MBax 
Parameter('kf13', transloc)
Parameter('kr13', transloc)
Rule('transloc_MBax_aBax', aBax(b=None) <> MBax(b=None), kf13, kr13)

# MBax + Bcl2 <-->  MBax:Bcl2  
Parameter('kf14', 1e-06/v)
Parameter('kr14', 1e-03)
Rule('bind_MBax_Bcl2', MBax(b=None) + Bcl2(b=None) <> MBax(b=1) * Bcl2(b=1), kf14, kr14)

# MBax + MBax <-->  Bax2
Parameter('kf15', 1e-06/v*2)
Parameter('kr15', 1e-03)
Rule('bind_MBax_MBax', MBax(b=None) + MBax(b=None) <> Bax2(b=None), kf15, kr15)

# Bax2 + Bcl2 <-->  Bax2:Bcl2  
Parameter('kf16', 1e-06/v)
Parameter('kr16', 1e-03)
Rule('bind_Bax2_Bcl2', Bax2(b=None) + Bcl2(b=None) <> Bax2(b=1) * Bcl2(b=1), kf16, kr16)

# Bax2 + Bax2 <-->  Bax4
Parameter('kf17', 1e-06/v*2)
Parameter('kr17', 1e-03)
Rule('bind_Bax2_Bax2', Bax2(b=None) + Bax2(b=None) <> Bax4(b=None), kf17, kr17)

# Bax4 + Bcl2 <-->  Bax4:Bcl2  
Parameter('kf18', 1e-06/v)
Parameter('kr18', 1e-03)
Rule('bind_Bax4_Bcl2', Bax4(b=None) + Bcl2(b=None) <> Bax4(b=1) * Bcl2(b=1), kf18, kr18)

# Bax4 + Mito <-->  Bax4:Mito -->  AMito  
Parameter('kf19', 1e-06/v)
Parameter('kr19', 1e-03)
Parameter('kc19', 1e+00)
Rule('bind_Bax4_Mito', Bax4(b=None) + Mito(b=None) <> Bax4(b=1) * Mito(b=1), kf19, kr19)
Rule('produce_AMito', Bax4(b=1) * Mito(b=1) >> AMito(b=None), kc19)

# AMito + mCytoC <-->  AMito:mCytoC --> AMito + ACytoC  
Parameter('kf20', 2e-06/v)
Parameter('kr20', 1e-03)
Parameter('kc20', 1e+01)
Rule('bind_AMito_mCytoC', AMito(b=None) + mCytoC(b=None) <> AMito(b=1) * mCytoC(b=1), kf20, kr20)
Rule('cat_AMito_ACytoC', AMito(b=1) * mCytoC(b=1) >> AMito(b=None) + ACytoC(b=None), kc20)

# AMito + mSmac <-->  AMito:mSmac --> AMito + ASmac  
Parameter('kf21', 2e-06/v)
Parameter('kr21', 1e-03)
Parameter('kc21', 1e+01)
Rule('bind_AMito_mSmac', AMito(b=None) + mSmac(b=None) <> AMito(b=1) * mSmac(b=1), kf21, kr21)
Rule('cat_AMito_ASmac', AMito(b=1) * mSmac(b=1) >> AMito(b=None) + ASmac(b=None), kc21)

# ACytoC <-->  cCytoC
Parameter('kf22', transloc)
Parameter('kr22', transloc)
Rule('transloc_cCytoC_ACytoC', ACytoC(b=None) <> cCytoC(b=None), kf22, kr22)

# Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
Parameter('kf23', 5e-07)
Parameter('kr23', 1e-03)
Parameter('kc23', 1e+00)
Rule('bind_Apaf_cCytoC', Apaf(b=None) + cCytoC(b=None) <> Apaf(b=1) * cCytoC(b=1), kf23, kr23)
Rule('cat_aApaf_cCytoC', Apaf(b=1) * cCytoC(b=1) >> aApaf(b=None) + cCytoC(b=None), kc23)

# aApaf + pC9 <-->  Apop
Parameter('kf24', 5e-08)
Parameter('kr24', 1e-03)
Rule('bind_aApaf_pC9', aApaf(b=None) + pC9(b=None) <> Apop(b=None), kf24, kr24)

# Apop + pC3 <-->  Apop:pC3 --> Apop + C3
Parameter('kf25', 5e-09)
Parameter('kr25', 1e-03)
Parameter('kc25', 1e+00)
Rule('bind_Apop_pC3', Apop(b=None) + pC3(b=None) <> Apop(b=1) * pC3(b=1), kf25, kr25)
Rule('cat_Apop_C3', Apop(b=1) * pC3(b=1) >> Apop(b=None) + C3(b=None), kc25)

# ASmac <-->  cSmac
Parameter('kf26', transloc)
Parameter('kr26', transloc)
Rule('transloc_cSmac_ASmac', ASmac(b=None) <> cSmac(b=None), kf26, kr26)

# Apop + XIAP <-->  Apop:XIAP  
Parameter('kf27', 2e-06)
Parameter('kr27', 1e-03)
Rule('bind_Apop_XIAP', Apop(b=None) + XIAP(b=None) <> Apop(b=1) * XIAP(b=1), kf27, kr27)

# cSmac + XIAP <-->  cSmac:XIAP  
Parameter('kf28', 7e-06)
Parameter('kr28', 1e-03)
Rule('bind_cSmac_XIAP', cSmac(b=None) + XIAP(b=None) <> cSmac(b=1) * XIAP(b=1), kf28, kr28)



# Fig 4B
observe('Bid',   Bid(b=None))
observe('PARP',  PARP(b=None))
observe('mSmac', mSmac(b=None))
# this is what I originally thought 4B was actually plotting
observe('tBid',  tBid())
observe('CPARP', CPARP())
observe('cSmac', cSmac())



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
        initial(m(sites), ic_param)














####


if __name__ == '__main__':
    from pysb.generator.bng import BngGenerator
    gen = BngGenerator(model)
    print gen.get_content()
    print ""
    print "begin actions"
    print "  generate_network({overwrite=>1});"
    print "  simulate_ode({t_end=>21600,n_steps=>360});" # 6 hours, 1-minute steps
    print "end actions"
