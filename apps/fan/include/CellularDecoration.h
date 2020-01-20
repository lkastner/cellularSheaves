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
#ifndef TROPICAL_DECORATION
#define TROPICAL_DECORATION

#include "polymake/client.h"
#include "polymake/graph/Decoration.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/BasicLatticeTypes.h"
#include "polymake/graph/lattice_builder.h"
#include "polymake/graph/LatticePermutation.h"
#include "polymake/fan/hasse_diagram.h"
#include "polymake/fan/CompactificationData.h"

namespace polymake { namespace fan{

   struct CellularDecoration : public GenericStruct<CellularDecoration> {
     DeclSTRUCT( DeclFIELD(face, Set<Int>)
                 DeclFIELD(rank,Int)
                 DeclFIELD(realisation, Set<Int>) );

     CellularDecoration() {}
     CellularDecoration(const Set<Int>& f, Int r, const Set<Int>& re) : face(f), rank(r), realisation(re) {}
   };
   
   class CellularDecorator {
      private:
         const Matrix<Rational>& vertices;
         const Map<Int, Set<Int>>& int2vertices;

         Set<Int> realisation(const Set<Int> face) const {
            Set<Int> result;
            for(const auto& e:face){
               result += int2vertices[e];
            }
            return result;
         }
      
      public:
         typedef CellularDecoration DecorationType;
         CellularDecorator(const CompactificationData& cd): vertices(cd.vertices), int2vertices(cd.int2vertices){}

         CellularDecoration compute_initial_decoration(const Set<Int>& face) const {
            vertices.rows();
            return CellularDecoration(face, 0, realisation(face));
         }

         CellularDecoration compute_decoration(const Set<Int>& face, const CellularDecoration& bd) const {
            return CellularDecoration(face, bd.rank+1, realisation(face));
         }

         CellularDecoration compute_artificial_decoration(const NodeMap<Directed, CellularDecoration>& decor, const std::list<Int>& max_faces) const {
            Set<Int> D;
            D += -1;
            return CellularDecoration(D, -1, D);
         }
         
   };


} // namespace fan
} // namespace polymake
#endif
