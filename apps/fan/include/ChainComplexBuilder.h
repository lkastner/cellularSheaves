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
            Map<Set<int>, Matrix<Rational>> bases;
            cosheaf.give("BLOCKS") >> blocksEM;
            cosheaf.give("BASES") >> bases;
            pc.give("ORIENTATIONS") >> orientationsEM;
            G = hd.graph();
            nodeDimsNM = NodeMap<Directed, int>(G);
            for(const auto& node:nodes(G)){
               nodeDimsNM[node] = bases[decoration[node].face].rows();
            }
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
                     result.minor(rowRanges[target], colRanges[source]) = orientationsEM[*edge] * blocksEM[*edge];
                  }
               }
            }
            return result;
         }
   };

   template<typename SelectorType>
   Array<Matrix<Rational>> build_chain_complex_from_hasse(perl::Object pc, perl::Object cosheaf, const SelectorType& selector, bool cochain){
      ChainComplexBuilder<SelectorType> SD(pc, cosheaf, selector);
      int dim = pc.give("FAN_DIM");
      Array<Matrix<Rational>> result(dim-1);
      for(int i=1; i<dim; i++){
         Matrix<Rational> B = SD.assemble_ith_matrix(i);
         if(cochain){
            result[(dim-1) - i] = T(B);
         } else {
            result[i-1] = B;
         }
      }
      return result;
   }
} // end namespace fan
} // end namespace polymake

#endif
