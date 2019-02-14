# GOMO
Generalized Operator Modelling of the Ocean (GOMO) is a practical ocean model based on OpenArray which is a simple operator library for the decoupling of ocean modelling and parallel computing.

The fundamental equations of GOMO are derived from POM. GOMO features a bottom-following, free-surface, staggered Arakawa C grid. To effectively evolve the rapid surface fluctuations, GOMO uses the mode-splitting algorithm to address the fast propagating surface gravity waves and slow propagating internal waves in barotropic (external) and baroclinic (internal) modes, respectively. The details of the continuous governing equations, the corresponding operator expression form and the descriptions of all the variables used in GOMO are listed in the docs folder.

To run GOMO, OpenArray should be compiled and installed at first, then compile GOMO through 'make' and run GOMO.
