// This is an automatically generated file, any manual changes done here will be overwritten!
// To change the generated C++ code, edit the script support/generate_cpperl_modules.pl
// To change the number of instances per source file, specify "chunksize":N in ../auto-choose_basis.cpperl

#define POLYMAKE_CPPERL_FILE "auto-choose_basis"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/client.h"
#include "polymake/fan/linalg_tools.h"
#include "polymake/linalg.h"

namespace polymake { namespace fan { namespace {
FunctionCallerStart4perl {
enum class choose_basis;
};
FunctionCaller4perl(choose_basis, free_t);
FunctionTemplateInstance4perl(0, choose_basis, free_t, choose_basis:T1.X, perl::Returns::normal, 1, (Rational, perl::Canned<const pm::RepeatedRow<pm::SameElementVector<pm::Rational const&> >&>));
FunctionTemplateInstance4perl(1, choose_basis, free_t, choose_basis:T1.X, perl::Returns::normal, 1, (Rational, perl::Canned<const Matrix<Rational>&>));
} } }
