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
#ifndef CHAIN_COMPLEX_BUILDER_H
#define CHAIN_COMPLEX_BUILDER_H

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


namespace polymake { namespace fan{
   
   template<typename NodeSelector>
   class ChainComplexBuilder {
      private: 
         Lattice<BasicDecoration, lattice::Nonsequential> hd;
         Graph<Directed> G;
         EdgeMap<Directed, Matrix<Rational>> blocksEM;
         EdgeMap<Directed, int> orientationsEM;
         NodeMap<Directed, int> nodeDimsNM;
         const NodeSelector& nodeSelector;
      public:
         ChainComplexBuilder(perl::Object pc, perl::Object cosheaf, const NodeSelector& nose):
            nodeSelector(nose)
         {
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

         bool node_is_valid(int node) const {
            return nodeSelector.isValid(hd.decoration()[node].face);
         }

         Matrix<Rational> assemble_ith_matrix(int i) const {
            int nrows=0, ncols=0;
            NodeMap<Directed, Set<int>> colRanges(G), rowRanges(G);
            for(const auto& node:hd.nodes_of_rank(i+1)){
               if(node_is_valid(node)){
                  colRanges[node] = pm::range(ncols, ncols+nodeDimsNM[node]-1);
                  // cout << node << ": " << nodeDimsNM[node] << " - " << colRanges[node] << endl;
                  ncols += nodeDimsNM[node];
               }
            }
            // cout << "cols done." << endl;
            for(const auto& node:hd.nodes_of_rank(i)){
               if(node_is_valid(node)){
                  rowRanges[node] = pm::range(nrows, nrows+nodeDimsNM[node]-1);
                  // cout << node << ": " << nodeDimsNM[node] << " - " << rowRanges[node] << endl;
                  nrows += nodeDimsNM[node];
               }
            }
            Matrix<Rational> result(nrows,ncols);
            int target;
            for(const auto& source: hd.nodes_of_rank(i+1)){
               for(auto edge = entire(G.in_edges(source)); !edge.at_end(); ++edge){
                  target = edge.from_node();
                  if(node_is_valid(source) && node_is_valid(target)){
                     // cout << target << " - " << source << endl;
                     // cout << rowRanges[target] << " - " << colRanges[source] << endl;
                     // cout << "Insert block: " << endl << blocksEM[*edge] << endl;
                     result.minor(rowRanges[target], colRanges[source]) = orientationsEM[*edge] * blocksEM[*edge];
                     // result.minor(rowRanges[target], colRanges[source]) = blocksEM[*edge];
                  }
               }
            }
            return result;
         }
   };

   template<typename SelectorType>
   Array<Matrix<Rational>> build_chain_complex_from_hasse(perl::Object pc, perl::Object cosheaf, const SelectorType& selector){
      ChainComplexBuilder<SelectorType> SD(pc, cosheaf, selector);
      int dim = pc.give("FAN_DIM");
      Array<Matrix<Rational>> result(dim-1);
      for(int i=1; i<dim; i++){
         Matrix<Rational> B = SD.assemble_ith_matrix(i);
         // if(cochain){
         //    result[(dim-1) - (i-1)] = B;
         // } else {
            result[i-1] = B;
         //}
      }
      return result;
   }
} // end namespace fan
} // end namespace polymake

#endif
