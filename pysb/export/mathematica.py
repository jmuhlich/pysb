"""
Module containing a class for converting a PySB model to a set of ordinary
differential equations for integration or analysis in Mathematica.

For information on how to use the model exporters, see the documentation
for :py:mod:`pysb.export`.

Output for the Robertson example model
======================================

The Mathematica code produced will follow the form as given below for
``pysb.examples.robertson``::

    (*
    A simple three-species chemical kinetics system known as "Robertson's
    example", as presented in:

    H. H. Robertson, The solution of a set of reaction rate equations, in Numerical
    Analysis: An Introduction, J. Walsh, ed., Academic Press, 1966, pp. 178-182.

    Mathematica model definition file for model robertson.
    Generated by pysb.export.mathematica.MathematicaExporter.

    Run with (for example):
    tmax = 10
    soln = NDSolve[Join[odes, initconds], slist, {t, 0, tmax}]
    Plot[s0[t] /. soln, {t, 0, tmax}, PlotRange -> All]
    *)

    (* Parameters *)
    k1 = 0.040000000000000001;
    k2 = 30000000;
    k3 = 10000;
    A0 = 1;
    B0 = 0;
    C0 = 0;

    (* List of Species *)
    (* s0[t] = A() *)
    (* s1[t] = B() *)
    (* s2[t] = C() *)

    (* ODEs *)
    odes = {
    s0'[t] == -k1*s0[t] + k3*s1[t]*s2[t],
    s1'[t] == k1*s0[t] - k2*s1[t]^2 - k3*s1[t]*s2[t],
    s2'[t] == k2*s1[t]^2
    }

    (* Initial Conditions *)
    initconds = {
    s0[0] == A0,
    s1[0] == B0,
    s2[0] == C0
    }

    (* List of Variables (e.g., as an argument to NDSolve) *)
    solvelist = {
    s0[t],
    s1[t],
    s2[t]
    }

    (* Run the simulation -- example *)
    tmax = 100
    soln = NDSolve[Join[odes, initconds], solvelist, {t, 0, tmax}]

    (* Observables *)
    Atotal = (s0[t] * 1) /. soln
    Btotal = (s1[t] * 1) /. soln
    Ctotal = (s2[t] * 1) /. soln

The output consists of a block of commands that define the ODEs, parameters,
species and other variables for the model, along with a set of descriptive
comments. The sections are as follows:

* The header comments identify the model and show an example of how to
  integrate the ODEs in Mathematica.
* The parameters block defines the numerical values of the named parameters.
* The list of species gives the mapping between the indexed species (``s0``,
  ``s1``, ``s2``) and their representation in PySB (``A()``, ``B()``, ``C()``).
* The ODEs block defines the set of ordinary differential equations and assigns
  the set of equations to the variable ``odes``.
* The initial conditions block defines the initial values for each species and
  assigns the set of conditions to the variable ``initconds``.
* The "list of variables" block enumerates all of the species in the model
  (``s0[t]``, ``s1[t]``, ``s2[t]``) and assigns them to the variable
  ``solvelist``; this list can be passed to the Mathematica command ``NDSolve``
  to indicate the variables to be solved for.
* This is followed by an example of how to call ``NDSolve`` to integrate the
  equations.
* Finally, the observables block enumerates the observables in the model,
  expressing each one as a linear combination of the appropriate species in
  the model. The interpolating functions returned by ``NDSolve`` are substituted
  in from the solution variable ``soln``, allowing the observables to be
  plotted.

Note that Mathematica does not permit underscores in variable names, so
any underscores used in PySB variables will be removed (e.g., ``A_total`` will
be converted to ``Atotal``).
"""

import pysb
import pysb.bng
import sympy
import re
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
from pysb.export import Exporter, ExpressionsNotSupported, \
    CompartmentsNotSupported


class MathematicaExporter(Exporter):
    """A class for returning the ODEs for a given PySB model for use in
    Mathematica.

    Inherits from :py:class:`pysb.export.Exporter`, which implements
    basic functionality for all exporters.
    """

    def export(self):
        """Generate the corresponding Mathematica ODEs for the PySB model
        associated with the exporter.

        Returns
        -------
        string
            String containing the Mathematica code for the model's ODEs.
        """
        if self.model.expressions:
            raise ExpressionsNotSupported()
        if self.model.compartments:
            raise CompartmentsNotSupported()

        output = StringIO()
        pysb.bng.generate_equations(self.model)

        # Add docstring if there is one
        if self.docstring:
            output.write('(*\n' + self.docstring + '\n')
        else:
            output.write("(*\n")

        # Header comment
        output.write("Mathematica model definition file for ")
        output.write("model " + self.model.name + ".\n")
        output.write("Generated by " \
                     "pysb.export.mathematica.MathematicaExporter.\n")
        output.write("\n")
        output.write("Run with (for example):\n")
        output.write("tmax = 10\n")
        output.write("soln = NDSolve[Join[odes, initconds], slist, " \
                     "{t, 0, tmax}]\n")
        output.write("Plot[s0[t] /. soln, {t, 0, tmax}, PlotRange -> All]\n")
        output.write("*)\n\n")

        # PARAMETERS
        # Note that in Mathematica, underscores are not allowed in variable
        # names, so we simply strip them out here
        params_str = ''
        for i, p in enumerate(self.model.parameters):
            # Remove underscores
            pname = p.name.replace('_', '')

            # Convert parameter values to scientific notation
            # If the parameter is 0, don't take the log!
            if p.value == 0:
                params_str += '%s = %g;\n' %  (pname, p.value)
            # Otherwise, take the log (base 10) and format accordingly
            else:
                val_str = '%.17g' % p.value
                if 'e' in val_str:
                    (mantissa, exponent) = val_str.split('e')
                    params_str += '%s = %s * 10^%s;\n' % \
                            (pname, mantissa, exponent)
                else:
                    params_str += '%s = %s;\n' %  (pname, val_str)

        ## ODEs ###
        odes_str = 'odes = {\n'
        # Concatenate the equations
        odes_str += ',\n'.join(['s%d == %s' %
                                (i, sympy.ccode(self.model.odes[i]))
                                for i in range(len(self.model.odes))])
        # Replace, e.g., s0 with s[0]
        odes_str = re.sub(r's(\d+)', lambda m: 's%s[t]' % (int(m.group(1))),
                          odes_str)
        # Add the derivative symbol ' to the left hand sides
        odes_str = re.sub(r's(\d+)\[t\] ==', r"s\1'[t] ==", odes_str)
        # Correct the exponentiation syntax
        odes_str = re.sub(r'pow\(([^,]+), ([^)]+)\)', r'\1^\2', odes_str)
        odes_str += '\n}'
        #c_code = odes_str
        # Eliminate underscores from parameter names in equations
        for i, p in enumerate(self.model.parameters):
            odes_str = re.sub(r'\b(%s)\b' % p.name, p.name.replace('_', ''),
                              odes_str)

        ## INITIAL CONDITIONS
        ic_values = ['0'] * len(self.model.odes)
        for i, ic in enumerate(self.model.initials):
            idx = self.model.get_species_index(ic.pattern)
            ic_values[idx] = ic.value.name.replace('_', '')

        init_conds_str = 'initconds = {\n'
        init_conds_str += ',\n'.join(['s%s[0] == %s' % (i, val)
                                     for i, val in enumerate(ic_values)])
        init_conds_str += '\n}'

        ## SOLVE LIST
        solvelist_str = 'solvelist = {\n'
        solvelist_str += ',\n'.join(['s%s[t]' % (i)
                                    for i in range(0, len(self.model.odes))])
        solvelist_str += '\n}'

        ## OBSERVABLES
        observables_str = ''
        for obs in self.model.observables:
            # Remove underscores
            observables_str += obs.name.replace('_', '') + ' = '
            #groups = self.model.observable_groups[obs_name]
            observables_str += ' + '.join(['(s%s[t] * %d)' % (s, c)
                               for s, c in zip(obs.species, obs.coefficients)])
            observables_str += ' /. soln\n' 

        # Add comments identifying the species
        species_str = '\n'.join(['(* s%d[t] = %s *)' % (i, s) for i, s in
                            enumerate(self.model.species)])

        output.write('(* Parameters *)\n')
        output.write(params_str + "\n")
        output.write('(* List of Species *)\n')
        output.write(species_str + "\n\n")
        output.write('(* ODEs *)\n')
        output.write(odes_str + "\n\n")
        output.write('(* Initial Conditions *)\n')
        output.write(init_conds_str + "\n\n")
        output.write('(* List of Variables (e.g., as an argument to NDSolve) ' \
                     '*)\n')
        output.write(solvelist_str + '\n\n')
        output.write('(* Run the simulation -- example *)\n')
        output.write('tmax = 100\n')
        output.write('soln = NDSolve[Join[odes, initconds], ' \
                     'solvelist, {t, 0, tmax}]\n\n')
        output.write('(* Observables *)\n')
        output.write(observables_str + '\n')

        return output.getvalue()

