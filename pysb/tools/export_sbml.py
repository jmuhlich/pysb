#!/usr/bin/env python

# FIXME this should use libsbml if available

import pysb
import pysb.bng
import sympy
import re
import sys
import os
from StringIO import StringIO

def run(model):
    output = StringIO()
    pysb.bng.generate_equations(model)

    output.write(
        """<?xml version="1.0" encoding="UTF-8"?>
  <sbml xmlns="http://www.sbml.org/sbml/level2/version4" xmlns:celldesigner="http://www.sbml.org/2001/ns/celldesigner" level="2" version="4">
    <model>
""")
    output.write(
        """<annotation>
  <celldesigner:extension>
  <celldesigner:modelVersion>4.0</celldesigner:modelVersion>
  <celldesigner:listOfSpeciesAliases>
""")
    for i in range(len(model.species)):
        output.write(
            """    <celldesigner:speciesAlias id="sa%d" species="s%d">
      <celldesigner:bounds x="0" y="0" w="80.0" h="40.0"/>
      <celldesigner:view state="usual"/>
    </celldesigner:speciesAlias>
""" % (i, i))
    output.write(
        """  </celldesigner:listOfSpeciesAliases>
  <celldesigner:listOfProteins>
""")
    for i, cp in enumerate(model.species):
        name = str(cp).replace('% ', '._br_')  # CellDesigner does something weird with % in names
        output.write('    <celldesigner:protein id="pr%d" name="%s" type="GENERIC"/>\n' % (i, name))
    output.write("""  </celldesigner:listOfProteins>
  </celldesigner:extension>
</annotation>
""")
    output.write(
        """        <listOfCompartments>
            <compartment id="default" name="default" spatialDimensions="0"/>
        </listOfCompartments>
""")

    ics = [[s, 0] for s in model.species]  # complexpattern, initial value
    for cp, ic_param in model.initial_conditions:
        ics[model.get_species_index(cp)][1] = ic_param.value
    output.write("        <listOfSpecies>\n")
    for i, (cp, value) in enumerate(ics):
        name = str(cp).replace('%', '.')  # CellDesigner does something weird with % in names
        output.write('            <species id="s%d" name="%s" compartment="default" initialAmount="%e">\n' % (i, name, value));
        output.write(
            """                <annotation>
                    <celldesigner:extension>
                      <celldesigner:speciesIdentity>
                        <celldesigner:class>PROTEIN</celldesigner:class>
                        <celldesigner:proteinReference>pr%d</celldesigner:proteinReference>
                      </celldesigner:speciesIdentity>
                    </celldesigner:extension>
                </annotation>
            </species>
""" % i)

    output.write("        </listOfSpecies>\n")

    output.write("        <listOfParameters>\n")
    for i, param in enumerate(model.parameters_rules()):
        output.write('            <parameter id="%s" name="%s" value="%e"/>\n' % (param.name, param.name, param.value));
    output.write("        </listOfParameters>\n")

    output.write("        <listOfReactions>\n")
    for i, reaction in enumerate(model.reactions_bidirectional):
        reversible = str(reaction['reversible']).lower()
        output.write('            <reaction id="r%d" name="r%d" reversible="%s">\n' % (i, i, reversible));
        output.write(
            """      <annotation>
        <celldesigner:extension>
          <celldesigner:reactionType>HETERODIMER_ASSOCIATION</celldesigner:reactionType>
          <celldesigner:baseReactants>
""")
        for i in reaction['reactants']:
            output.write('            <celldesigner:baseReactant species="s%d" alias="sa%d"/>\n' % (i, i))
        output.write("          </celldesigner:baseReactants>\n          <celldesigner:baseProducts>\n")
        for i in reaction['products']:
            output.write('            <celldesigner:baseProduct species="s%d" alias="sa%d"/>\n' % (i, i))
        output.write(
            """          </celldesigner:baseProducts>
          <celldesigner:editPoints>0,0</celldesigner:editPoints>
        </celldesigner:extension>
      </annotation>
""")

        output.write('                <listOfReactants>\n');
        for species in reaction['reactants']:
            output.write('                    <speciesReference species="s%d">\n' % species)
            output.write("""          <annotation>
            <celldesigner:extension>
              <celldesigner:alias>sa%d</celldesigner:alias>
            </celldesigner:extension>
          </annotation>
""" % species)
            output.write('                    </speciesReference>\n')
        output.write('                </listOfReactants>\n');
        output.write('                <listOfProducts>\n');
        for species in reaction['products']:
            output.write('                    <speciesReference species="s%d">\n' % species)
            output.write("""          <annotation>
            <celldesigner:extension>
              <celldesigner:alias>sa%d</celldesigner:alias>
            </celldesigner:extension>
          </annotation>
""" % species)
            output.write('                    </speciesReference>\n')
        output.write('                </listOfProducts>\n');
        formula = sympy.ccode(reaction['rate'])
        output.write('                <kineticLaw formula="%s"/>\n' % formula);
        output.write('            </reaction>\n');
    output.write("        </listOfReactions>\n")

    output.write("    </model>\n</sbml>\n")
    return output.getvalue()

if __name__ == '__main__':
    # sanity checks on filename
    if len(sys.argv) <= 1:
        raise Exception("You must specify the filename of a model script")
    model_filename = sys.argv[1]
    if not os.path.exists(model_filename):
        raise Exception("File '%s' doesn't exist" % model_filename)
    if not re.search(r'\.py$', model_filename):
        raise Exception("File '%s' is not a .py file" % model_filename)
    sys.path.insert(0, os.path.dirname(model_filename))
    model_name = re.sub(r'\.py$', '', os.path.basename(model_filename))
    # import it
    try:
        # FIXME if the model has the same name as some other "real" module
        # which we use, there will be trouble
        # (use the imp package and import as some safe name?)
        model_module = __import__(model_name)
    except StandardError as e:
        print "Error in model script:\n"
        raise
    # grab the 'model' variable from the module
    try:
        model = model_module.__dict__['model']
    except KeyError:
        raise Exception("File '%s' isn't a model file" % model_filename)
    print run(model)



