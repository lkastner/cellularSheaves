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


#include "polymake/client.h"
#include "polymake/GenericMatrix.h"
#include "polymake/Map.h"
#include "polymake/Set.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"
#include "polymake/fan/CompactificationData.h"
#include "polymake/fan/CellularClosure.h"
#include "polymake/fan/CellularDecoration.h"
#include "polymake/fan/SedentarityDecoration.h"
#include "polymake/fan/AugmentedHasseDiagram.h"
#include "polymake/fan/linalg_tools.h"
#include "polymake/ListMatrix.h"



namespace polymake { namespace fan{
 



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

   class SheafData {
      private: 
         Lattice<BasicDecoration, lattice::Nonsequential> hd;
         Graph<Directed> G;
         EdgeMap<Directed, Matrix<Rational>> blocksEM;
         EdgeMap<Directed, int> orientationsEM;
         NodeMap<Directed, int> nodeDimsNM;
      public:
         SheafData(perl::Object pc, perl::Object cosheaf){
            pc.give("HASSE_DIAGRAM") >> hd;
            const auto& decoration = hd.decoration();
            Map<Set<Set<int>>, Matrix<Rational>> blocks;
            Map<Set<Set<int>>, int> orientations;
            Map<Set<int>, Matrix<Rational>> bases;
            cosheaf.give("BLOCKS") >> blocks;
            cosheaf.give("BASES") >> bases;
            pc.give("ORIENTATIONS") >> orientations;
            G = hd.graph();
            blocksEM = EdgeMap<Directed, Matrix<Rational>>(G);
            orientationsEM = EdgeMap<Directed, int>(G);
            nodeDimsNM = NodeMap<Directed, int>(G);
            // cout << "Edges: " << edges(G) << endl;
            for(const auto& node:nodes(G)){
               nodeDimsNM[node] = bases[decoration[node].face].rows();
            }
            Set<int> source, target;
            Set<Set<int>> blocksKey;
            for(auto edge=entire(edges(G)); !edge.at_end(); ++edge){
               // cout << "E:" << *edge << ": " << decoration[edge.from_node()] << " - " << decoration[edge.to_node()] << endl;
               source = decoration[edge.from_node()].face;
               target = decoration[edge.to_node()].face;
               blocksKey.clear();
               blocksKey += source;
               blocksKey += target;
               // cout << blocks[blocksKey] << endl;
               blocksEM[*edge] = blocks[blocksKey];
               if(orientations.exists(blocksKey)){
                  // cout << "Has orientation: " << orientations[blocksKey] << endl;
                  orientationsEM[*edge] = orientations[blocksKey];
               }
               if(nodeDimsNM[edge.from_node()] == blocks[blocksKey].rows()){
                  // cout << "rows";
               }
               if(nodeDimsNM[edge.from_node()] == blocks[blocksKey].cols()){
                  // cout << "cols";
               }
            }
            // cout << nodeDimsNM << endl;
         }

         Matrix<Rational> assemble_ith_matrix(int i) const {
            int nrows=0, ncols=0;
            NodeMap<Directed, Set<int>> colRanges(G), rowRanges(G);
            for(const auto& node:hd.nodes_of_rank(i+1)){
               colRanges[node] = pm::range(ncols, ncols+nodeDimsNM[node]-1);
               // cout << node << ": " << nodeDimsNM[node] << " - " << colRanges[node] << endl;
               ncols += nodeDimsNM[node];
            }
            // cout << "cols done." << endl;
            for(const auto& node:hd.nodes_of_rank(i)){
               rowRanges[node] = pm::range(nrows, nrows+nodeDimsNM[node]-1);
               // cout << node << ": " << nodeDimsNM[node] << " - " << rowRanges[node] << endl;
               nrows += nodeDimsNM[node];
            }
            Matrix<Rational> result(nrows,ncols);
            int target;
            for(const auto& source: hd.nodes_of_rank(i+1)){
               for(auto edge = entire(G.in_edges(source)); !edge.at_end(); ++edge){
                  target = edge.from_node();
                  // cout << target << " - " << source << endl;
                  // cout << rowRanges[target] << " - " << colRanges[source] << endl;
                  // cout << "Insert block: " << endl << blocksEM[*edge] << endl;
                  result.minor(rowRanges[target], colRanges[source]) = orientationsEM[*edge] * blocksEM[*edge];
                  // result.minor(rowRanges[target], colRanges[source]) = blocksEM[*edge];
               }
            }
            return result;
         }
   };

   void check_complex(perl::Object pc, perl::Object cosheaf, bool cochain){
      if(!cochain){ return; }
      SheafData SD(pc, cosheaf);
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
         Matrix<Rational> A = assemble_matrix_cpp(sigmas, taus, blocks, orientations);
         Matrix<Rational> B = SD.assemble_ith_matrix(i);
         cout << A.rows() << " " << A.cols() << endl;
         cout << B.rows() << " " << B.cols() << endl;
         cout << "Check: " << (A==B) << endl;

         // cout << assemble_matrix_cpp(sigmas, taus, blocks, orientations) << endl;
      }
      
   }

   Function4perl(&tropcomp, "tropcomp( $ )");
   
   Function4perl(&check_complex, "check_complex( $ , $ , $ )");

} // namespace fan
} // namespace polymake
