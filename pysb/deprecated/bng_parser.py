import sys
import argparse
import warnings
from ply import lex, yacc;
import re
import pysb.core
import pysb.export


reserved_list = [
    'begin',
    'end',
    'model',
    'parameters',
    'molecule',
    'types',
    'observables',
    'seed',
    'species',
    'reaction',
    'rules',
    'reactions',
    'groups',
    'substanceUnits',
    ]
reserved = dict([r, r.upper()] for r in reserved_list)

tokens = [
    'ID',
    
    'FLOAT',
    'INTEGER',
    'STRING',

    'COMMA',
    'PLUS',
    'ASTERISK',
    'EXCLAMATION',
    'AT',
    'QUESTION',
    'PERIOD',
    'COLON',
    'DOLLAR',
    'SLASH',
    'EQUALS',
    'ARROW_REVERSIBLE',
    'ARROW_IRREVERSIBLE',
    'LPAREN',
    'RPAREN',
    'LBRACE',
    'RBRACE',
    'ARROW_ARGS',
    'SEMI',
    'SITE_STATE',
    'NL',
    ] + reserved.values()

t_COMMA = r','
t_PLUS = r'\+'
t_ASTERISK = r'\*'
t_EXCLAMATION = r'!'
t_AT = r'@'
t_QUESTION = r'\?'
t_PERIOD = r'\.'
t_COLON = r':'
t_DOLLAR = r'\$'
t_SLASH = r'/'
t_EQUALS = r'='
t_ARROW_IRREVERSIBLE = r'->'
t_ARROW_REVERSIBLE = r'<->'
t_LPAREN = r'\('
t_RPAREN = r'\)'
t_LBRACE = r'{'
t_RBRACE = r'}'
t_ARROW_ARGS = r'=>'
t_SEMI = r';'

def t_STRING(t):
    r'"[^"]*"'
    t.value = t.value[1:-1]
    return t

# Define a rule so we can track line numbers
def t_NL(t):
    r'(\n|\r|\r\n)'
    t.lexer.lineno += 1
    return t

def t_FLOAT(t):
    r'[+-]?(\d*\.\d+([eE][+-]?\d+)?|\d+[eE][+-]?\d+)'
    try:
        t.value = float(t.value)    
    except ValueError:
        print "Line %d: Number '%s' has some kind of problem (ValueError)!" % (t.lineno,t.value)
        t.value = float("nan")
    return t

def t_INTEGER(t):
    r'\d+'
    try:
        t.value = int(t.value)    
    except ValueError:
        print "Line %d: Number '%s' has some kind of problem (ValueError)!" % (t.lineno,t.value)
        t.value = 0
    return t

def t_ID(t):
    r'[a-zA-Z_][a-zA-Z_0-9]*'
    t.type = reserved.get(t.value,'ID')  # check for reserved words
    return t

def t_SITE_STATE(t):
    r'~[A-Za-z0-9_]+'
    t.value = t.value[1:]
    return t

# Match and ignore comments (# to end of line)
def t_comment(t):
    r'\#[^\n\r]*'

def t_line_cont(t):
    r'\\\s*(\n|\r|\r\n)'

# A string containing ignored characters (spaces and tabs)
t_ignore  = ' \t'

# Error handling rule
def t_error(t):
    print "Illegal character '%s' on line %d" % (t.value[0], t.lineno)
    t.lexer.skip(1)


def list_helper(p):
    if len(p) == 1:
        p[0] = []
    if len(p) == 2:
        p[0] = [p[1]]
    elif len(p) == 3:
        p[0] = p[1] + [p[2]]
    p[0] = [v for v in p[0] if v != None] # filter out Nones

def comma_list_helper(p):
    if len(p) == 1:
        p[0] = []
    if len(p) == 2:
        p[0] = [p[1]]
    elif len(p) == 4:
        p[0] = p[1] + [p[3]]


def p_model(p):
    '''model : block_list'''
    p[0] = p.parser.pysb_model

def p_block_list(p):
    '''block_list : block_list block
                  | block'''

def p_block(p):
    '''block : substance_units
             | parameter_block
             | molecule_type_block
             | species_block
             | reaction_rules_block
             | reactions_block
             | observables_block
             | block_empty'''

def p_block_empty(p):
    '''block_empty : BEGIN MODEL eol
                   | END MODEL eol
                   | eol'''

def p_substance_units(p):
    '''substance_units : SUBSTANCEUNITS LPAREN STRING RPAREN eol'''

def p_parameter_block(p):
    'parameter_block : BEGIN PARAMETERS eol parameter_st_list END PARAMETERS eol'

def p_molecule_type_block(p):
    'molecule_type_block : BEGIN MOLECULE TYPES eol molecule_type_st_list END MOLECULE TYPES eol'

def p_species_block(p):
    '''species_block : BEGIN SPECIES eol END SPECIES eol
                     | BEGIN SEED SPECIES eol END SEED SPECIES eol'''

def p_reaction_rules_block(p):
    'reaction_rules_block : BEGIN REACTION RULES eol reaction_rule_st_list END REACTION RULES eol'

def p_reactions_block(p):
    'reactions_block : BEGIN REACTIONS eol END REACTIONS eol'

def p_observables_block(p):
    'observables_block : BEGIN OBSERVABLES eol END OBSERVABLES eol'


def p_parameter_st_list(p):
    '''parameter_st_list : parameter_st_list parameter_st
                         | parameter_st'''
    list_helper(p)

def p_parameter_st(p):
    '''parameter_st : INTEGER parameter_dec
                    | parameter_dec'''
    if len(p) == 3:
        p[0] = p[2]
    elif len(p) == 2:
        p[0] = p[1]

def p_parameter_dec(p):
    '''parameter_dec : ID number eol'''
    parameter = pysb.core.Parameter(p[1], p[2], _export=False)
    p.parser.pysb_model.add_component(parameter)
    #print '====>', parameter


def p_molecule_type_st_list(p):
    '''molecule_type_st_list : molecule_type_st_list molecule_type_st
                             | molecule_type_st'''
    list_helper(p)

def p_molecule_type_st(p):
    '''molecule_type_st : INTEGER molecule_type_dec
                        | molecule_type_dec'''
    if len(p) == 3:
        p[0] = p[2]
    else:
        p[0] = p[1]

def p_molecule_type_dec(p):
    '''molecule_type_dec : ID LPAREN site_def_list RPAREN eol'''
    sites = [s[0] for s in p[3]]
    site_states = dict((k, v) for k, v in p[3] if v is not None)
    monomer = pysb.core.Monomer(p[1], sites, site_states, _export=False)
    p.parser.pysb_model.add_component(monomer)
    #print '====>', monomer

def p_site_def_list(p):
    '''site_def_list : site_def_list COMMA site_def
                     | site_def
                     |'''
    comma_list_helper(p)

def p_site_def(p):
    '''site_def : ID site_state_list
                | ID'''
    if len(p) == 2:
        p[0] = (p[1], None)
    elif len(p) == 3:
        p[0] = (p[1], p[2])

def p_site_state_list(p):
    '''site_state_list : site_state_list SITE_STATE
                       | SITE_STATE'''
    list_helper(p)


def p_reaction_rule_st_list(p):
    '''reaction_rule_st_list : reaction_rule_st_list reaction_rule_st
                             | reaction_rule_st'''
    list_helper(p)

def p_reaction_rule_st(p):
    '''reaction_rule_st : ID COLON rule_expression rule_rates eol
                        | rule_expression rule_rates eol'''
    if len(p) == 4:
        expr, rates = p[1], p[2]
        num_rules = len(p.parser.pysb_model.rules)
        name = 'Rule' + str(num_rules + 1)
    elif len(p) == 6:
        name, expr, rates = p[1], p[3], p[4]
    rule = pysb.core.Rule(name, expr, *rates, _export=False)
    p.parser.pysb_model.add_component(rule)
    #print '====>', rule

def p_rule_expression(p):
    '''rule_expression : rule_expression_reversible
                       | rule_expression_irreversible'''
    p[0] = p[1]

def p_rule_expression_reversible(p):
    '''rule_expression_reversible : reaction_pattern ARROW_REVERSIBLE reaction_pattern'''
    p[0] = pysb.core.RuleExpression(p[1], p[3], True)

def p_rule_expression_irreversible(p):
    '''rule_expression_irreversible : reaction_pattern ARROW_IRREVERSIBLE reaction_pattern'''
    p[0] = pysb.core.RuleExpression(p[1], p[3], False)

def p_reaction_pattern(p):
    '''reaction_pattern : reaction_pattern PLUS complex_pattern
                        | complex_pattern'''
    if len(p) == 2:
        p[0] = pysb.core.ReactionPattern([p[1]])
    elif len(p) == 4:
        p[0] = p[1] + p[3]

def p_complex_pattern(p):
    '''complex_pattern : complex_pattern PERIOD monomer_pattern
                       | monomer_pattern'''
    if len(p) == 2:
        p[0] = pysb.core.ComplexPattern([p[1]], None)
    elif len(p) == 4:
        p[0] = p[1] % p[3]

def p_monomer_pattern(p):
    '''monomer_pattern : ID LPAREN site_condition_list RPAREN'''
    try:
        monomer = p.parser.pysb_model.monomers[p[1]]
    except KeyError:
        import ipdb; ipdb.set_trace()
        raise SyntaxError("Unknown molecule type: %s" % p[1])
    p[0] = monomer(dict(p[3]))
    

def p_site_condition_list(p):
    '''site_condition_list : site_condition_list COMMA site_condition
                           | site_condition
                           |'''
    comma_list_helper(p)

def p_site_condition(p):
    '''site_condition : ID SITE_STATE bond
                      | ID bond
                      | ID SITE_STATE
                      | ID'''
    if len(p) == 2:
        p[0] = (p[1], None)
    elif len(p) == 3:
        p[0] = (p[1], p[2])
    elif len(p) == 4:
        p[0] = (p[1], (p[2], p[3]))

def p_bond(p):
    '''bond : EXCLAMATION bond_state'''
    p[0] = p[2]

def p_bond_state(p):
    '''bond_state : INTEGER
                  | PLUS
                  | QUESTION'''
    if p[1] == '+':
        p[0] = pysb.core.ANY
    elif p[1] == '?':
        p[0] = pysb.core.WILD
    else:
        p[0] = p[1]


def p_rule_rates(p):
    '''rule_rates : ID COMMA ID
                  | ID'''
    if len(p) == 2:
        names = [p[1]]
    elif len(p) == 4:
        names = [p[1], p[3]]
    p[0] = [p.parser.pysb_model.parameters[name] for name in names]

def p_number(p):
    '''number : FLOAT
              | INTEGER'''
    p[0] = p[1]

# This was easier than making the NL regexp match multiple newlines. That would
# have entailed properly counting line numbers in the face of CR vs. LF vs. CRLF
# line endings.
def p_eol(p):
    '''eol : eol NL
           | NL'''


# Error rule for syntax errors
def p_error(p):
    print "Syntax error on line %d:" % p.lineno
    print p


precedence = (
    ('left', 'PLUS'),
)


block_names = ['parameter', 'molecule_type', 'observable', 'species',
               'reaction_rule', 'reaction', 'group']

class Block(object):
    """AST-type node to hold the contents of a begin/end block."""
    def __init__(self, name, contents):
        if name not in block_names:
            raise InvalidBlockNameError(name)
        self.name = name
        self.contents = contents
    def __repr__(self):
        return 'Block(%s, <%d items>)' % (repr(self.name), len(self.contents))

class InvalidBlockNameError(ValueError):
    """Block name is invalid."""


# Build the lexer and parser
lexer = lex.lex()
parser = yacc.yacc(write_tables=0)
parser.pysb_model = pysb.core.Model(name='model', _export=False)


def parse(*args, **kwargs):
    yacc.parse(*args, **kwargs)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("-l", "--lexer", help="test lexer", action='store_true')
    ap.add_argument("-p", "--parser", help="test parser", action='store_true')
    ap.add_argument("-d", "--debug", help="debug output", action='store_true')
    ap.add_argument("bngl_file", nargs="?")
    args = ap.parse_args()
    if args.bngl_file is None:
        from pysb.examples.robertson import model
        content = pysb.export.export(model, 'bng_net')
    else:
        with open(args.bngl_file) as f:
            content = f.read()
    logger = None
    if args.debug:
        logger = yacc.PlyLogger(sys.stdout)
    if args.lexer:
        lexer.input(content)
        for tok in lexer:
            print tok
    elif args.parser:
        model = parser.parse(content, debug=logger)
        content = pysb.export.export(model, 'bngl')
        print "\n======================\n"
        print content
    else:
        print "ERROR: Must specify -l or -p"
        ap.print_help()
