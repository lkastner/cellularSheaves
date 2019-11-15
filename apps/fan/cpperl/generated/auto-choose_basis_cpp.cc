// This is an automatically generated file, any manual changes done here will be overwritten!
// To change the generated C++ code, edit the script support/generate_cpperl_modules.pl
// To change the number of instances per source file, specify "chunksize":N in ../auto-choose_basis_cpp.cpperl

#define POLYMAKE_CPPERL_FILE "auto-choose_basis_cpp"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Set.h"
#include "polymake/client.h"
#include "polymake/fan/linalg_tools.h"

namespace polymake { namespace fan { namespace {
FunctionCallerStart4perl {
enum class choose_basis_cpp;
};
FunctionCaller4perl(choose_basis_cpp, free_t);
FunctionTemplateInstance4perl(0, choose_basis_cpp, free_t, choose_basis_cpp:T1.X, perl::Returns::normal, 1, (Rational, perl::Canned<const Matrix<Rational>&>));
FunctionTemplateInstance4perl(1, choose_basis_cpp, free_t, choose_basis_cpp:T1.X, perl::Returns::normal, 1, (Rational, perl::Canned<const pm::BlockMatrix<mlist<pm::MatrixMinor<pm::Matrix<pm::Rational> const&, pm::Set<int, pm::operations::cmp> const&, pm::all_selector const&> const&, pm::Matrix<pm::Rational> const&>, std::integral_constant<bool, true> >&>));
} } }
