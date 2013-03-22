import sys
import argparse
from ply import lex, yacc;
import pysb
import pysb.core


reserved_list = [
    'begin',
    'end',
    'parameters',
    'molecule',
    'types',
    'species',
    'reaction',
    'rules',
    'reactions',
    'observables',
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
    'TILDE',
    'EXCLAMATION',
    'QUESTION',
    'PERIOD',
    'COLON',
    'ARROW_REVERSIBLE',
    'ARROW_IRREVERSIBLE',
    'LPAREN',
    'RPAREN',
    'NEWLINE',
    ] + reserved.values()

t_COMMA = r','
t_PLUS = r'\+'
t_TILDE = r'~'
t_EXCLAMATION = r'!'
t_QUESTION = r'\?'
t_PERIOD = r'\.'
t_COLON = r':'
t_ARROW_IRREVERSIBLE = r'-->'
t_ARROW_REVERSIBLE = r'<->'
t_LPAREN = r'\('
t_RPAREN = r'\)'

def t_STRING(t):
    r'"[^"]*"'
    t.value = t.value[1:-1]
    return t

# Define a rule so we can track line numbers
def t_NEWLINE(t):
    r'\n+'
    t.lexer.lineno += len(t.value)
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

# Match and ignore comments (# to end of line)
def t_comment(t):
    r'\#[^\n]*'

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
    'model : block_list'
    model = Model(name='model', _export=False)
    p[0] = model

def p_block_list(p):
    '''block_list : block_list block
                  | block'''
    list_helper(p)

def p_block(p):
    '''block : substance_units
             | parameter_block
             | molecule_type_block
             | species_block
             | reaction_rules_block
             | observables_block'''
    p[0] = p[1]

def p_block_empty(p):
    '''block : NEWLINE'''

def p_parameter_block(p):
    'parameter_block : BEGIN PARAMETERS NEWLINE parameter_st_list END PARAMETERS NEWLINE'
    p[0] = pysb.core.ComponentSet(p[4])

def p_molecule_type_block(p):
    'molecule_type_block : BEGIN MOLECULE TYPES NEWLINE molecule_type_st_list END MOLECULE TYPES NEWLINE'
    p[0] = p[5]

def p_molecule_type_st_list(p):
    '''molecule_type_st_list : molecule_type_st_list molecule_type_st
                             | molecule_type_st
                             | NEWLINE'''
    list_helper(p)

def p_molecule_type_st(p):
    '''molecule_type_st : INTEGER molecule_type_dec
                        | molecule_type_dec'''
    if len(p) == 3:
        p[0] = p[2]
    else:
        p[0] = p[1]

def p_molecule_type_dec(p):
    'molecule_type_dec : ID LPAREN site_def_list RPAREN'
    p[0] = pysb.Monomer(p[1], _export=False)

def p_site_def_list(p):
    '''site_def_list : site_def_list COMMA site_def
                     | site_def
                     |'''
    comma_list_helper(p)

def p_site_def(p):
    '''site_def : ID'''
    p[0] = p[1]

def p_species_block(p):
    'species_block : BEGIN SPECIES NEWLINE END SPECIES NEWLINE'
    p[0] = p[2]

def p_reaction_rules_block(p):
    'reaction_rules_block : BEGIN REACTION RULES NEWLINE END REACTION RULES NEWLINE'
    p[0] = p[2]

def p_observables_block(p):
    'observables_block : BEGIN OBSERVABLES NEWLINE END OBSERVABLES NEWLINE'
    p[0] = p[2]

def p_parameter_st_list(p):
    '''parameter_st_list : parameter_st_list parameter_st
                         | parameter_st
                         | NEWLINE'''
    list_helper(p)

def p_parameter_st(p):
    '''parameter_st : INTEGER parameter_dec
                    | parameter_dec'''
    if len(p) == 3:
        p[0] = p[2]
    else:
        p[0] = p[1]

def p_parameter_dec(p):
    '''parameter_dec : ID number NEWLINE'''
    p[0] = pysb.Parameter(p[1], p[2], _export=False)

def p_number(p):
    '''number : FLOAT
              | INTEGER'''
    p[0] = p[1]

def p_substance_units(p):
    '''substance_units : SUBSTANCEUNITS LPAREN STRING RPAREN NEWLINE'''

# def p_statement_list(p):
#     '''statement_list : statement_list statement'''
#     p[0] = p[1] + [p[2]]
#     #print "statement_list:", p[0]

# def p_statement_list_trivial(p):
#     '''statement_list : statement'''
#     p[0] = [p[1]]
#     #print "statement_list_trivial:", p[0]

# def p_statement_empty(p):
#     'statement : NEWLINE'
#     #print "statement_empty:", p[0]

# def p_statement(p):
#     'statement : rule NEWLINE'
#     p[0] = p[1]
#     #print "statement:", p[0]

# def p_rule(p):
#     '''rule : irr_rule
#             | rev_rule'''
#     p[0] = p[1]
#     #print "rule:", p[0]

# def p_irr_rule(p):
#     'irr_rule : expression IRRARROW expression LPAREN FLOAT RPAREN'
#     #print "irr_rule"
#     p[0] = RuleIrreversible(reactants=p[1], products=p[3], rate=p[5])

# def p_rev_rule(p):
#     'rev_rule : expression REVARROW expression LPAREN FLOAT COMMA FLOAT RPAREN'
#     #print "rev_rule"
#     p[0] = RuleReversible(reactants=p[1], products=p[3], rates=[p[5], p[7]])

# def p_expression_plus(p):
#     'expression : expression PLUS expression'
#     p[0] = p[1] + p[3]
#     #print "expression_plus:", p[0]

# def p_expression_species(p):
#     'expression : SPECIES'
#     #print "expression_species:", p[1]
#     p[0] = [Species(name=p[1])]

# Error rule for syntax errors
def p_error(p):
    print "Syntax error in input:"
    print p

precedence = (
    ('left', 'PLUS'),
)

# Build the lexer and parser
lexer = lex.lex()
parser = yacc.yacc(write_tables=0)


def parse(*args, **kwargs):
    yacc.parse(*args, **kwargs)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("-l", "--lexer", help="debug lexer", action='store_true')
    ap.add_argument("-p", "--parser", help="debug parser", action='store_true')
    ap.add_argument("bngl_file", nargs="?")
    args = ap.parse_args()
    if args.bngl_file is None:
        from pysb.examples.robertson import model
        from pysb.export import export
        content = export(model, 'bng_net')
    else:
        with open(args.bngl_file) as f:
            content = f.read()
    if args.lexer:
        lexer.input(content)
        for tok in lexer:
            print tok
    elif args.parser:
        parser.parse(content, debug=yacc.PlyLogger(sys.stdout))
    else:
        print "ERROR: Must specify -l or -p"
        ap.print_help()
