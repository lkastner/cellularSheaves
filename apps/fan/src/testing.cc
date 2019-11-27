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


#include "polymake/fan/ChainComplexBuilder.h"


namespace polymake { namespace fan{

   void check_complex(perl::Object pc, perl::Object cosheaf, bool cochain){
      /*
      if(!cochain){ return; }
      TrivialSelector ts;
      BoundedSelector bs(pc);
      // ChainComplexBuilder<TrivialSelector> SD(pc, cosheaf, ts);
      ChainComplexBuilder<BoundedSelector> SD(pc, cosheaf, bs);
      Lattice<BasicDecoration, lattice::Nonsequential> hd = pc.give("HASSE_DIAGRAM");
      Map<Set<Set<int>>, Matrix<Rational>> blocks;
      Map<Set<Set<int>>, int> orientations;
      Map<Set<int>, Matrix<Rational>> bases;
      cosheaf.give("BLOCKS") >> blocks;
      pc.give("ORIENTATIONS") >> orientations;
      // cout << "Blocks: " << blocks << endl;
      // cout << "Orientations: " << orientations << endl;
      // cout << "Decoration" << endl << hd.decoration() << endl;
      // cout << hd.nodes_of_rank(1) << endl;
      int dim = pc.give("FAN_DIM");
      const auto& decoration = hd.decoration();
      for(int i=1; i<dim; i++){
         // cout << i << ": " << hd.nodes_of_rank(i) << endl;
         Array<Set<int>> sigmas(hd.nodes_of_rank(i+1).size());
         int count = 0;
         for(const auto s:hd.nodes_of_rank(i+1)){
            sigmas[count] = decoration[s].face;
            count++;
         }
         Array<Set<int>> taus(hd.nodes_of_rank(i).size());
         count = 0;
         for(const auto s:hd.nodes_of_rank(i)){
            taus[count] = decoration[s].face;
            count++;
         }
         // cout << "Sigmas: " << sigmas << endl;
         // cout << "Taus: " << taus << endl;
         // Matrix<Rational> A = assemble_matrix_cpp(sigmas, taus, blocks, orientations);
         // Matrix<Rational> B = SD.assemble_ith_matrix(i);
         // cout << B << endl;
         // cout << A.rows() << " " << A.cols() << endl;
         // cout << B.rows() << " " << B.cols() << endl;
         // cout << "Check: " << (A==B) << endl;

         // cout << assemble_matrix_cpp(sigmas, taus, blocks, orientations) << endl;
      }
      */
      
   }
   
   Function4perl(&check_complex, "check_complex( $ , $ , $ )");

} // namespace fan
} // namespace polymake
