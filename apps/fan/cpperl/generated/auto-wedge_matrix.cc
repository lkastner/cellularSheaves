// This is an automatically generated file, any manual changes done here will be overwritten!
// To change the generated C++ code, edit the script support/generate_cpperl_modules.pl
// To change the number of instances per source file, specify "chunksize":N in ../auto-wedge_matrix.cpperl

#define POLYMAKE_CPPERL_FILE "auto-wedge_matrix"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/client.h"
#include "polymake/fan/linalg_tools.h"

namespace polymake { namespace fan { namespace {
FunctionCallerStart4perl {
enum class wedge_matrix;
};
FunctionCaller4perl(wedge_matrix, free_t);
FunctionTemplateInstance4perl(0, wedge_matrix, free_t, wedge_matrix:T1.X.x, perl::Returns::normal, 1, (Rational, perl::Canned<const Matrix<Rational>&>, void));
FunctionTemplateInstance4perl(1, wedge_matrix, free_t, wedge_matrix:T1.X.x, perl::Returns::normal, 1, (Rational, perl::Canned<const pm::Transposed<pm::Matrix<pm::Rational> >&>, void));
} } }
