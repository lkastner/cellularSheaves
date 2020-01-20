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
#ifndef TROPICAL_CLOSURE
#define TROPICAL_CLOSURE

#include "polymake/client.h"
#include "polymake/graph/Decoration.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/BasicLatticeTypes.h"
#include "polymake/graph/lattice_builder.h"
#include "polymake/graph/LatticePermutation.h"
#include "polymake/fan/hasse_diagram.h"
#include "polymake/fan/CellularDecoration.h"
#include "polymake/fan/CompactificationData.h"

namespace polymake { namespace fan{
	
   using graph::Lattice;
   using namespace graph::lattice;
   using namespace fan::lattice;
   
   struct IteratorWrap {
      FacetList groundSet;
      FacetList::const_iterator state;

      IteratorWrap(FacetList&& gs): groundSet(gs){
         state = groundSet.begin();
      }

      Set<Int> operator*(){
         return *state;
      }

      void operator++(){
         state++;
      }

      bool at_end(){
         return state == groundSet.end();
      }
   };

   template<typename DecorationType>
   class CellularClosureOperator {
      private:
         FaceMap<> face_index_map;
         const Map<Int, Set<Int>>& int2vertices;
         const Map<Set<Int>, Int>& vertices2int;
         const Int& nVertices;
         const Set<Int>& farVertices;
         ComplexPrimalClosure<> oldClosureOperator;

         IncidenceMatrix<> construct_old_closure_operator(perl::Object pc) {
            RestrictedIncidenceMatrix<> building_matrix;
            const Array<IncidenceMatrix<> >& maximal_vifs = pc.give("MAXIMAL_CONES_INCIDENCES");
            bool is_pure, is_complete = false;
            pc.give("PURE") >> is_pure;
            TopologicalType tt(is_pure, is_complete);
            const IncidenceMatrix<>& maximal_cones = pc.give("MAXIMAL_CONES");
            const Int n_vertices = maximal_cones.cols();
            FacetList non_redundant_facets(n_vertices);
            for (auto mvf : maximal_vifs) {
               for (auto fct = entire(rows(mvf)); !fct.at_end(); ++fct)
                  non_redundant_facets.replaceMax(*fct);
            }

            building_matrix /= maximal_cones;
            for (auto nrf = entire(non_redundant_facets); !nrf.at_end(); ++nrf)
               building_matrix /= *nrf;

            return IncidenceMatrix<>(std::move(building_matrix));
         }

      public:
         typedef Set<Int> ClosureData;

         CellularClosureOperator(const CompactificationData& cd, perl::Object pc):
            int2vertices(cd.int2vertices), vertices2int(cd.vertices2int), nVertices(cd.nVertices), farVertices(cd.farVertices), oldClosureOperator(construct_old_closure_operator(pc)){
            }
         
         Set<Int> old_closure(const Set<Int>& a) const {
            Set<Int> result(range(0,nVertices));
            bool contained = false;
            for(const auto& pair: vertices2int){
               if(incl(a, pair.first)<=0){
                  contained = true;
                  result = result*pair.first;
               }
            }
            if(!contained){
               result.clear();
               result += -1;
            }
            return result;
         }

         Set<Int> closure(const Set<Int> a) const {
            Set<Int> originalRealisation;
            for(const auto i:a){
               originalRealisation += int2vertices[i];
            }
            Set<Int> originalClosure = old_closure(originalRealisation);
            Set<Int> commonRays = originalRealisation * farVertices;
            for(const auto i : a){
               commonRays = commonRays * int2vertices[i];
            }
            Set<Int> result;
            for(const auto& v:vertices2int){
               if(incl(commonRays, v.first)<=0 && incl(v.first, originalClosure)<=0){
                  result += v.second;
               }
            }
            return result;
         }
         
         Set<Int> closure_of_empty_set(){
            Set<Int> empty;
            return empty;
         }

         FaceIndexingData get_indexing_data(const ClosureData& data)
         {
            Int& fi = face_index_map[data];
            return FaceIndexingData(fi, fi == -1, fi == -2);
         }

         Set<Int> compute_closure_data(const DecorationType& bd) const {
            return bd.face;
         }

         IteratorWrap get_closure_iterator(const Set<Int>& face) const {
            Set<Int> all = pm::range(0,int2vertices.size()-1);
            Set<Int> toadd = all-face;
            FacetList result;
            for(auto i:toadd){
               result.insertMin(closure(face+i));
            }
            return IteratorWrap(std::move(result));
         }
   };

} // namespace fan
} // namespace polymake
#endif
