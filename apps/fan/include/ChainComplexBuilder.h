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
   
   template<typename HasseDiagram, typename NodeSelector, typename MatrixType>
   class ChainComplexBuilder {
      private: 
         const HasseDiagram& hd;
         const Graph<Directed>& G;
         const EdgeMap<Directed, MatrixType>& orientedBlocks;
         const NodeMap<Directed, Int>& nodeDimsNM;
         const NodeSelector& nodeSelector;
      public:
         ChainComplexBuilder(const HasseDiagram& hd_in,
            const Graph<Directed>& G_in,
            const EdgeMap<Directed, MatrixType>& ob_in,
            const NodeMap<Directed, Int>& nd_in,
            const NodeSelector& nose):
            hd(hd_in), G(G_in), orientedBlocks(ob_in), nodeDimsNM(nd_in), nodeSelector(nose)
         {}

         bool node_is_valid(Int node) const {
            return nodeSelector.isValid(hd.decoration()[node].face);
         }

         MatrixType assemble_ith_matrix(Int i) const {
            Int nrows=0, ncols=0;
            NodeMap<Directed, Set<Int>> colRanges(G), rowRanges(G);
            for(const auto& node:hd.nodes_of_rank(i)){
               if(node_is_valid(node)){
                  colRanges[node] = pm::range(ncols, ncols+nodeDimsNM[node]-1);
                  ncols += nodeDimsNM[node];
               }
            }
            for(const auto& node:hd.nodes_of_rank(i+1)){
               if(node_is_valid(node)){
                  rowRanges[node] = pm::range(nrows, nrows+nodeDimsNM[node]-1);
                  nrows += nodeDimsNM[node];
               }
            }
            MatrixType result(nrows,ncols);
            Int target;
            for(const auto& source: hd.nodes_of_rank(i+1)){
               for(auto edge = entire(G.in_edges(source)); !edge.at_end(); ++edge){
                  target = edge.from_node();
                  if(node_is_valid(source) && node_is_valid(target)){
                     result.minor(rowRanges[source], colRanges[target]) = orientedBlocks[*edge];
#if POLYMAKE_DEBUG
                     MatrixType insert(orientedBlocks[*edge]);
                     if(rowRanges[source].size() != insert.rows() || colRanges[target].size() != insert.cols()){
                        throw std::runtime_error("Matrix dimension does not agree with minor dimension.");
                     }
#endif
                  }
               }
            }
            return result;
         }
   };
   
   template<typename SelectorType, typename HasseDiagramType, typename MatrixType>
   topaz::ChainComplex<MatrixType> build_chain_complex_from_hasse(const HasseDiagramType& hd, const EdgeMap<Directed, Int> orientations, const EdgeMap<Directed, MatrixType>& blocks_in, const SelectorType& selector, bool cochain){
      const Graph<Directed>& G(hd.graph());
      NodeMap<Directed, Int> nodeDimsNM(G);
      EdgeMap<Directed, MatrixType> blocks(blocks_in);
      for(auto edge=entire(edges(G)); !edge.at_end(); ++edge){
         if(cochain){
            blocks[*edge] = T(blocks[*edge]);
         }
         nodeDimsNM[edge.from_node()] = blocks[*edge].cols();
         if(orientations[*edge] == -1){
            blocks[*edge] = -blocks[*edge];
         }
      }

      ChainComplexBuilder<HasseDiagramType, SelectorType, MatrixType> SD(hd, G, blocks, nodeDimsNM, selector);
      Int dim = hd.rank() - 1;
      Array<MatrixType> result(dim-1);
      for(Int i=1; i<dim; i++){
         result[i-1] = SD.assemble_ith_matrix(i);
      }
#if POLYMAKE_DEBUG
      return topaz::ChainComplex<MatrixType>(result, true);
#else
      return topaz::ChainComplex<MatrixType>(result);
#endif
   }

} // end namespace fan
} // end namespace polymake

#endif
