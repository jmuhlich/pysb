import re
from pprint import pprint


f = file('albeck_plosbio_2008_earm10.txt')
params = []
monomers = []
rules = []
monomers_seen = {}

mp_re = r'([A-Za-z0-9_:*]+)'
rp_re = r'%(mp_re)s *(?:\+ *%(mp_re)s)*' % {'mp_re': mp_re}
eq_re = r'(<?-->|==)'
#rule_re = r'%(rp_re)s *(?:(<?-->|==) *%(rp_re)s)+' % {'rp_re': rp_re}

param_names = []

for line in f:
    if 'rate constants' in line.lower():
        break
for line in f:
    if 'figure simulations' in line.lower():
        break
    elif line.startswith('k'):
        matches = re.finditer(r'([a-z_]+)\((\d+)\)=([0-9.E-]+|transloc)', line)
        param_names = []
        param_values = []
        for m in matches:
            (base, num, value) = m.groups()
            if base == 'k':
                base = 'kf'
            elif base == 'k_':
                base = 'kr'
            name = base + num
            try:
                value = '%.0e' % float(value)
            except ValueError:
                pass
            param_names.append(name)
            param_values.append(value)
            rules.append("Parameter('%s', %s)" % (name, value))
        rule_parts = re.split(eq_re, prev_rule_line.replace(' ', ''))
        species = rule_parts[0::2]
        eqs = rule_parts[1::2]
        param_pos = 0
        for ei in range(len(eqs)):
            eq = eqs[ei]
            reactants = species[ei].split('+')
            products = species[ei + 1].split('+')
            if eq == '-->':
                operator = '>>'
                eq_params = param_names[param_pos:param_pos+1]
                eq_param_values = param_values[param_pos:param_pos+1]
                param_pos += 1
            elif eq == '<-->':
                operator = '<>'
                eq_params = param_names[param_pos:param_pos+2]
                eq_param_values = param_values[param_pos:param_pos+2]
                param_pos += 2
            else:
                operator = None
            if len(reactants) == 2 and len(products) == 1:
                rule_name = 'bind_%s_%s' % tuple(reactants)
            elif len(reactants) == 1 and len(products) == 2:
                rule_name = 'cat_%s_%s' % tuple(products)
            elif len(reactants) == 1 and len(products) == 1 and eq_param_values[0] == 'transloc':
                rule_name = 'transloc_%s_%s' % (products[0], reactants[0])
            elif len(reactants) == 1 and len(products) == 1 and eq_param_values[0] <> 'transloc':
                rule_name = 'produce_%s' % (products[0])
            for s_list in reactants, products:
                for si, sp in enumerate(s_list):
                    if ':' in sp:
                        s_list[si] = ' * '.join(['%s(b=1)' % m for m in sp.split(':')])
                    else:
                        s_list[si] = '%s(b=None)' % sp
            r_str = ' + '.join(reactants)
            p_str = ' + '.join(products)
            rule = "Rule('%s', %s %s %s, %s)" % (rule_name, r_str, operator, p_str, ', '.join(eq_params))
            rules.append(rule)
            if ei == 0:
                rules.insert(-len(param_names) - 1, '\n#' + prev_rule_line)
        if param_pos <> len(param_names):
            print "ERROR: param list over- or underconsumed! pp:", param_pos, "len(pn):", len(param_names)
    elif line.startswith('% '):
        matches = re.finditer(mp_re, line)
        for m in matches:
            name = m.groups()[0]
            name = name.replace('*', '_act')
            if ':' in name:
                continue
            if name in monomers_seen:
                continue
            monomers_seen[name] = None
            monomers.append("Monomer('%s', ['b'])" % (name,))
        prev_rule_line = line.replace('%', '').replace('\n', '')

for m in monomers:
    print m
print
for r in rules:
    print r
