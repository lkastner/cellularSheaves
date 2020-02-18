// This is an automatically generated file, any manual changes done here will be overwritten!
// To change the generated C++ code, edit the script support/generate_cpperl_modules.pl

#define POLYMAKE_CPPERL_FILE "wrap-build_complex"
#include "polymake/perl/macros.h"
#include FindDefinitionSource4perl(build_complex.cc)
#include "polymake/GF2.h"
#include "polymake/Graph.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/fan/compactification.h"
#include "polymake/graph/Decoration.h"

namespace polymake { namespace fan { namespace {
FunctionCallerStart4perl {
enum class build_full_chain;
};
FunctionCaller4perl(build_full_chain, free_t);
FunctionTemplateInstance4perl(0, build_full_chain, free_t, build_full_chain:T3.B.X.X.x, perl::Returns::normal, 3, (fan::compactification::SedentarityDecoration, graph::lattice::Nonsequential, Rational, void, perl::Canned<const EdgeMap<Directed, Int>&>, perl::Canned<const EdgeMap<Directed, Matrix<Rational>>&>, void));
FunctionTemplateInstance4perl(1, build_full_chain, free_t, build_full_chain:T3.B.X.X.x, perl::Returns::normal, 3, (fan::compactification::SedentarityDecoration, graph::lattice::Nonsequential, GF2, void, perl::Canned<const EdgeMap<Directed, Int>&>, perl::Canned<const EdgeMap<Directed, Matrix<GF2>>&>, void));
} } }
