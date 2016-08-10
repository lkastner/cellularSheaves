/* Copyright (c) 1997-2014
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#ifndef POLYMAKE_WEDGE_MATRIX
#define POLYMAKE_WEDGE_MATRIX

#include "polymake/client.h"
#include "polymake/SparseVector.h"
#include "polymake/SparseMatrix.h"
#include "polymake/GenericMatrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/Map.h"
#include "polymake/Set.h"
#include "polymake/Integer.h"
#include "polymake/linalg.h"
#include "polymake/Array.h"
#include "polymake/PowerSet.h"
#include <list>

namespace polymake { namespace fan{
   

template <typename E> inline
Matrix<E>
wedge_matrix(const Matrix<E>& input, int n){
   int nrows, ncols, resultNrows, resultNcols, i, j;
   Matrix<E> result;
   nrows = input.rows();
   ncols = input.cols();
   resultNrows = (int) Integer::binom(nrows, n);
   resultNcols = (int) Integer::binom(ncols, n);
   if((resultNrows == 0) || (resultNcols == 0)){
      return zero_matrix<E>(resultNrows, resultNcols);
   }
   result = Matrix<E>(resultNrows, resultNcols);
   i=0;
   for (auto selectedRow=entire(all_subsets_of_k(sequence(0,nrows),n)); !selectedRow.at_end(); ++selectedRow) {
      j=0;
      cout << *selectedRow << endl;
      for (auto selectedCol=entire(all_subsets_of_k(sequence(0,ncols),n)); !selectedCol.at_end(); ++selectedCol) {
         result(i,j) = det(input.minor(*selectedRow, *selectedCol));
         j++;
      }
      i++;
   };
   return result;
}

template <typename E> inline
Matrix<E>
choose_basis(const Matrix<E>& input){
	int desired = rank(input), current = 0, test;
	Matrix<E> result(0, input.cols());
	for(auto rowit = entire(rows(input)); !rowit.at_end(); ++rowit){
		if(desired == rank(result)){
			return result;
		}
		test = rank(result / *rowit);
		if(test > current){
			current = test;
			result /= *rowit;
		}
	}
   return input;
}


} // namespace fan
} // namespace polymake
#endif
