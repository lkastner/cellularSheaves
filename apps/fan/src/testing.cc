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
         NonFarSelector(perl::Object pc){
            pc.give("FAR_VERTICES") >> farFace;
         }

         bool isValid(const Set<int>& face) const{
            return !(face - farFace).empty();
         }
   };
   
   class BoundedSelector {
      private:
         Set<int> farFace;
      public:
         BoundedSelector(perl::Object pc){
            pc.give("FAR_VERTICES") >> farFace;
         }

         bool isValid(const Set<int>& face) const{
            return (face * farFace).empty();
         }
   };


   // perl::Object
   void tropcomp(perl::Object pc){
      CompactificationData cd(pc);
      CellularDecorator decorator(cd);
      AugmentedHasseDiagram<CellularDecorator> AHD(cd, decorator, pc);
      AHD.print();

      Matrix<Rational> O(ones_matrix<Rational>(10,10));
      cout << O;

      pm::perl::PropertyValue vert(pc.give("VERTICES"));
      cout << "V:" << endl;
      
      //std::cout << vert << endl;
      std::ostringstream buffer;
      // vert.put(std::forward<std::ostringstream>(buffer));
      // vert.parse(buffer);
      // auto wrapped_buffer = wrap(buffer);
      // // wrapped_buffer << polymake::legible_typename(typeid(obj)) << pm::endl;

      // wrapped_buffer << O;
      std::cout << buffer.str();

      // Use a different decorator:
      cout << endl << "-------------------------------" << endl;
      cout << "Using SedentarityDecorator: " << endl;
      SedentarityDecorator sd(cd);
      AugmentedHasseDiagram<SedentarityDecorator> AHDSD(cd, sd, pc);
      AHDSD.print();
    
      // const Lattice<CellularDecoration, lattice::Nonsequential>& HD(AHD.get_HasseDiagram());
      // cout << "Graph: " << endl;
      // cout << HD.graph() << endl;
      // cout << "Node decoration: " << endl;
      // cout << HD.decoration() << endl;
      // return HD.makeObject();
      
     //  template <typename Decoration, typename ClosureOperator, typename CrossCut, typename Decorator, bool dual, typename SeqType = lattice::Nonsequential >
     //  Lattice<Decoration, SeqType> compute_lattice_from_closure(
     //        ClosureOperator cl,
     //        const CrossCut& cut,
     //        const Decorator& decorator,
     //        bool wants_artificial_top_node,
     //        bool_constant<dual> built_dually,
     //        Lattice<Decoration, SeqType> lattice = Lattice<Decoration>(),
     //        Set<int> queuing_nodes = Set<int>())
      ListMatrix<Vector<Rational>> A;
      A /= zero_vector<Rational>(5);
      A /= zero_vector<Rational>(5);
      A /= ones_vector<Rational>(5);
      cout << *(--(rows(A).end())) << endl;
   }


   void check_complex(perl::Object pc, perl::Object cosheaf, bool cochain){
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
      cosheaf.give("BASES") >> bases;
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
      
   }
   
   Array<Matrix<Rational>> build_nonfar_chain(perl::Object pc, perl::Object cosheaf, bool cochain){
      // if(!cochain){ return; }
      NonFarSelector bs(pc);
      // ChainComplexBuilder<TrivialSelector> SD(pc, cosheaf, ts);
      Array<Matrix<Rational>> result(build_chain_complex_from_hasse(pc, cosheaf, bs, cochain));
      // cout << "Returning" << endl;
      return result;
   }

   
   Array<Matrix<Rational>> build_bounded_chain(perl::Object pc, perl::Object cosheaf, bool cochain){
      // if(!cochain){ return; }
      BoundedSelector bs(pc);
      // ChainComplexBuilder<TrivialSelector> SD(pc, cosheaf, ts);
      Array<Matrix<Rational>> result(build_chain_complex_from_hasse(pc, cosheaf, bs, cochain));
      // cout << "Returning" << endl;
      return result;
   }

   Function4perl(&tropcomp, "tropcomp( $ )");
   
   Function4perl(&check_complex, "check_complex( $ , $ , $ )");
   
   Function4perl(&build_bounded_chain, "build_bounded_chain( $ , $ , $ )");
   
   Function4perl(&build_nonfar_chain, "build_nonfar_chain( $ , $ , $ )");

} // namespace fan
} // namespace polymake
