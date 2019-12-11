/* Copyright (c) 2016-2019
   Lars Kastner (TU Berlin)
   Kristin Shaw (University of Oslo)
   Anna-Lena Winz (FU Berlin)

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
      // cout << *selectedRow << endl;
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
   if(desired == 0){
      // If the matrix is empty
      return result;
   }
	for(const auto& rowit : rows(input)){
		if(desired == rank(result)){
			return result;
		}
		test = rank(result / rowit);
		if(test > current){
			current = test;
			result /= rowit;
		}
	}
   return result;
}


template <typename E> inline
Matrix<E>
build_matrix(const Matrix<E>& bigger, const Matrix<E>& smaller){
   // If one of the matrices is empty.
   if(bigger.rows() == 0){
      return zero_matrix<E>(smaller.rows(), 0);
   }
   if(smaller.rows() == 0){
      return zero_matrix<E>(0, bigger.rows());
   }
   Matrix<E> result(smaller.rows(), bigger.rows()+1),test,image;
   E lastVal;
   int j = 0;
   for(const auto& rowit : rows(smaller)){
      test = bigger / rowit;
      image = null_space(T(test));
      int i = 0;
      lastVal = image(i,image.cols()-1);
      while(lastVal == 0){
         i++;
         lastVal = image(i,image.cols()-1);
      }
      image = -1/lastVal * image;
      result.row(j) = image.row(i);
      j++;
   }
   return result.minor(All,~scalar2set(result.cols()-1));
}
 
   

} // namespace fan
} // namespace polymake
#endif
