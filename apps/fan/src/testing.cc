/* Copyright (c) 2016-2018
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
 
   typedef Lattice<BasicDecoration, lattice::Nonsequential> HasseDiagramType;

   class TrivialSelector {
      public:
         TrivialSelector(){}

         bool isValid(const Set<int>& face) const{
            return true;
         }
   };
   
   class NonFarSelector {
      private:
         Set<int> farFace;
      public:
         NonFarSelector(const Set<int>& ff) : farFace(ff) {}

         bool isValid(const Set<int>& face) const{
            return !(face - farFace).empty();
         }
   };
   
   class BoundedSelector {
      private:
         Set<int> farFace;
      public:
         BoundedSelector(const Set<int>& ff) : farFace(ff) {}

         bool isValid(const Set<int>& face) const{
            return (face * farFace).empty();
         }
   };




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
   
   Array<Matrix<Rational>> build_nonfar_chain(const HasseDiagramType& hd, const EdgeMap<Directed, int>& orientations, const Set<int>& ff, perl::Object cosheaf, bool cochain){
      NonFarSelector bs(ff);
      Array<Matrix<Rational>> result(build_chain_complex_from_hasse(hd, orientations, cosheaf, bs, cochain));
      return result;
   }

   
   Array<Matrix<Rational>> build_bounded_chain(const HasseDiagramType& hd, const EdgeMap<Directed, int>& orientations, const Set<int>& ff, perl::Object cosheaf, bool cochain){
      BoundedSelector bs(ff);
      Array<Matrix<Rational>> result(build_chain_complex_from_hasse(hd, orientations, cosheaf, bs, cochain));
      return result;
   }

   
   Function4perl(&check_complex, "check_complex( $ , $ , $ )");
   
   Function4perl(&build_bounded_chain, "build_bounded_chain( $ , $ , $ , $ , $ )");
   
   Function4perl(&build_nonfar_chain, "build_nonfar_chain( $ , $ , $ , $ , $ )");

} // namespace fan
} // namespace polymake
