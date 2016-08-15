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
wedge_matrix_cpp(const Matrix<E>& input, int n){
   int nrows, ncols, resultNrows, resultNcols, i, j;
   Matrix<E> result;
   nrows = input.rows();
   ncols = input.cols();
   // cout << "choose: " << ncols << " " << n<<endl;
   resultNrows = Integer::binom(nrows, n).to_int();
   resultNcols = Integer::binom(ncols, n).to_int();
   if((resultNrows == 0) || (resultNcols == 0)){
      return zero_matrix<E>(resultNrows, resultNcols);
   }
   result = Matrix<E>(resultNrows, resultNcols);
   i=0;
   for (Entire<Subsets_of_k<const sequence&> >::const_iterator selectedRow=entire(all_subsets_of_k(sequence(0,nrows),n)); !selectedRow.at_end(); ++selectedRow) {
      j=0;
      // cout << *selectedRow << endl;
      for (Entire<Subsets_of_k<const sequence&> >::const_iterator selectedCol=entire(all_subsets_of_k(sequence(0,ncols),n)); !selectedCol.at_end(); ++selectedCol) {
         result(i,j) = pm::det(input.minor(*selectedRow, *selectedCol));
         j++;
      }
      i++;
   };
   return result;
}

template <typename E> inline
Matrix<E>
choose_basis_cpp(const Matrix<E>& input){
	int desired = rank(input), current = 0, test;
	Matrix<E> result(0, input.cols());
	for(typename Entire<Rows<Matrix<E> > >::const_iterator rowit = entire(rows(input)); !rowit.at_end(); ++rowit){
		if(desired == rank(result)){
			return result;
		}
		test = rank(result / *rowit);
		if(test > current){
			current = test;
			result /= *rowit;
		}
	}
   return result;
}

template <typename E> inline
Matrix<E>
assemble_matrix_cpp(const Array<Set<int> > sigmas, 
               const Array<Set<int> > taus, 
               const Map<Set<Set<int> >, Matrix<E> > blocks, 
               const Map<Set<Set<int> >, int> orientations){
   int nrows = taus.size(), ncols = sigmas.size(), i, j, blockRow, blockCol;
   cout << "AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" << endl;
   Set<Set<int> > first_pair(sigmas[0]), currentPair;
   first_pair += Set<Set<int> >(taus[0]);
   cout << "FP: " << first_pair << endl;
   cout << blocks << endl;
   blockRow = blocks[first_pair].rows();
   blockCol = blocks[first_pair].cols();
   nrows *= blockRow;
   ncols *= blockCol;
   if((nrows == 0) || (ncols == 0)){
      return zero_matrix<E>(nrows, ncols);
   }
   Matrix<E> result(nrows, ncols);
   i = 0;
   for(typename Entire<Array<Set<int> > >::const_iterator tau = entire(taus); !tau.at_end(); ++tau){
      j = 0;
      for(typename Entire<Array<Set<int> > >::const_iterator sigma = entire(sigmas); !sigma.at_end(); ++sigma){
         currentPair = Set<Set<int> >(*sigma);
         currentPair += Set<Set<int> >(*tau);
         // cout << "Hello." << endl;
         // cout << result << endl;
         // cout << blocks[currentPair] << endl;
         // cout << i*blockRow << " " << (i+1)*blockRow << endl;
         // cout << j*blockCol << " " << (j+1)*blockCol << endl;
         // cout << sequence(3,5) << endl;
         result.minor(sequence(i*blockRow, blockRow), sequence(j*blockCol, blockCol)) = orientations[currentPair] * blocks[currentPair];
         j++;
      }
      i++;
   }
   return result;
   
}

template <typename E> inline
Matrix<E>
build_matrix_cpp(const Matrix<E>& bigger, const Matrix<E>& smaller){
  Matrix<E> result(0, bigger.cols()+1),test,image;
  E lastVal;
   for(typename Entire<Rows<Matrix<E> > >::const_iterator rowit = entire(rows(smaller)); !rowit.at_end(); ++rowit){
      test = bigger / *rowit;
      image = null_space(T(test));
      if (image.rows() != 1){
         
      }
      lastVal = image(0,image.cols()-1);
      image = -1/lastVal * image;
      result /= image;

   }
   return result.minor(All,~scalar2set(result.cols()-1));
}
 
   

} // namespace fan
} // namespace polymake
#endif
