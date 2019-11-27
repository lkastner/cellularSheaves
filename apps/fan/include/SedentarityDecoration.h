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
#ifndef SEDENTARITY_DECORATION
#define SEDENTARITY_DECORATION

#include "polymake/client.h"
#include "polymake/graph/Decoration.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/BasicLatticeTypes.h"
#include "polymake/graph/lattice_builder.h"
#include "polymake/graph/LatticePermutation.h"
#include "polymake/fan/hasse_diagram.h"
#include "polymake/fan/CompactificationData.h"

namespace polymake { namespace fan{

   struct SedentarityDecoration : public GenericStruct<SedentarityDecoration> {
     DeclSTRUCT( DeclFIELD(face, Set<int>)
                 DeclFIELD(rank,int)
                 DeclFIELD(realisation, Set<int>) 
                 DeclFIELD(sedentarity, Set<int>) );

     SedentarityDecoration() {}
     SedentarityDecoration(const Set<int>& f, int r, const Set<int>& re, const Set<int>& se) : face(f), rank(r), realisation(re), sedentarity(se) {}
   };
   
   class SedentarityDecorator {
      private:
         const Map<int, Set<int>>& int2vertices;
         const Set<int>& farVertices;

         Set<int> realisation(const Set<int> face) const {
            Set<int> result;
            for(const auto& e:face){
               result += int2vertices[e];
            }
            return result;
         }

         Set<int> sedentarity(const Set<int>& face) const {
            if(face.size() == 0){
               return Set<int>();
            }
            Set<int> result(farVertices);
            for(const auto& e:face){
               result *= int2vertices[e];
            }
            return result;
         }
      
      public:
         typedef SedentarityDecoration DecorationType;
         SedentarityDecorator(const CompactificationData& cd): int2vertices(cd.int2vertices), farVertices(cd.farVertices){}

         SedentarityDecoration compute_initial_decoration(const Set<int>& face) const {
            return SedentarityDecoration(face, 0, realisation(face), sedentarity(face));
         }

         SedentarityDecoration compute_decoration(const Set<int>& face, const SedentarityDecoration& bd) const {
            return SedentarityDecoration(face, bd.rank+1, realisation(face), sedentarity(face));
         }

         SedentarityDecoration compute_artificial_decoration(const NodeMap<Directed, SedentarityDecoration>& decor, const std::list<int>& max_faces) const {
            Set<int> D;
            Set<int> empty;
            D += -1;
            return SedentarityDecoration(D, -1, D, empty);
         }
         
   };


} // namespace fan
} // namespace polymake
#endif
