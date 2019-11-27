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
#include "polymake/topaz/ChainComplex.h"


namespace polymake { namespace fan{
   
   template<typename HasseDiagram, typename NodeSelector>
   class ChainComplexBuilder {
      private: 
         const HasseDiagram& hd;
         const Graph<Directed>& G;
         const EdgeMap<Directed, Matrix<Rational>>& orientedBlocks;
         const NodeMap<Directed, int>& nodeDimsNM;
         const NodeSelector& nodeSelector;
      public:
         ChainComplexBuilder(const HasseDiagram& hd_in,
            const Graph<Directed>& G_in,
            const EdgeMap<Directed, Matrix<Rational>>& ob_in,
            const NodeMap<Directed, int>& nd_in,
            const NodeSelector& nose):
            hd(hd_in), G(G_in), orientedBlocks(ob_in), nodeDimsNM(nd_in), nodeSelector(nose)
         {}

         bool node_is_valid(int node) const {
            return nodeSelector.isValid(hd.decoration()[node].face);
         }

         Matrix<Rational> assemble_ith_matrix(int i) const {
            int nrows=0, ncols=0;
            NodeMap<Directed, Set<int>> colRanges(G), rowRanges(G);
            for(const auto& node:hd.nodes_of_rank(i+1)){
               if(node_is_valid(node)){
                  colRanges[node] = pm::range(ncols, ncols+nodeDimsNM[node]-1);
                  ncols += nodeDimsNM[node];
               }
            }
            for(const auto& node:hd.nodes_of_rank(i)){
               if(node_is_valid(node)){
                  rowRanges[node] = pm::range(nrows, nrows+nodeDimsNM[node]-1);
                  nrows += nodeDimsNM[node];
               }
            }
            Matrix<Rational> result(nrows,ncols);
            int target;
            for(const auto& source: hd.nodes_of_rank(i+1)){
               for(auto edge = entire(G.in_edges(source)); !edge.at_end(); ++edge){
                  target = edge.from_node();
                  if(node_is_valid(source) && node_is_valid(target)){
                     // cout << "Prev: " << endl << result;
                     result.minor(rowRanges[target], colRanges[source]) = orientedBlocks[*edge];
                     Matrix<Rational> insert(orientedBlocks[*edge]);
                     if(rowRanges[target].size() != insert.rows() || colRanges[source].size() != insert.cols()){
                        throw std::runtime_error("Matrix dimension does not agree with minor dimension.");
                     }
                     // cout << "Inserting:" << endl << insert;
                     // cout << "After: " << endl << result;
                     // cout << "Rows: " << rowRanges[target] << " Cols: " << colRanges[source] << endl;
                     // cout << "iRows: " << insert.rows() << " iCols: " << insert.cols() << endl << endl;
                  }
               }
            }
            return result;
         }
   };
   
   template<typename SelectorType, typename HasseDiagramType>
   topaz::ChainComplex<Matrix<Rational>> build_chain_complex_from_hasse(const HasseDiagramType& hd, const EdgeMap<Directed, int> orientations, perl::Object cosheaf, const SelectorType& selector, bool cochain){
      const Graph<Directed>& G(hd.graph());
      EdgeMap<Directed, Matrix<Rational>> blocks;
      NodeMap<Directed, int> nodeDimsNM(G);
      cosheaf.give("BLOCKS") >> blocks;
      for(auto edge=entire(edges(G)); !edge.at_end(); ++edge){
         if(cochain){
            blocks[*edge] = T(blocks[*edge]);
         }
         nodeDimsNM[edge.from_node()] = blocks[*edge].rows();
         blocks[*edge] *= orientations[*edge];
         // if(!cochain){
         //    nodeDimsNM[edge.from_node()] = blocks[*edge].rows();
         // } else {
         //    nodeDimsNM[edge.from_node()] = blocks[*edge].cols();
         // }
         // cout << blocks[*edge] << endl;
         // blocks[*edge] = T(blocks[*edge]);
         // cout << blocks[*edge] << endl;
      }

      ChainComplexBuilder<HasseDiagramType, SelectorType> SD(hd, G, blocks, nodeDimsNM, selector);
      int dim = hd.rank() - 1;
      Array<Matrix<Rational>> result(dim-1);
      for(int i=1; i<dim; i++){
         Matrix<Rational> B = SD.assemble_ith_matrix(i);
         result[i-1] = T(B);
      }
      return topaz::ChainComplex<Matrix<Rational>>(result, true);
   }

} // end namespace fan
} // end namespace polymake

#endif
