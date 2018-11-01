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
#ifndef TROPICAL_DECORATION
#define TROPICAL_DECORATION

#include "polymake/client.h"
#include "polymake/graph/Decoration.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/BasicLatticeTypes.h"
#include "polymake/graph/lattice_builder.h"
#include "polymake/graph/LatticePermutation.h"
#include "polymake/fan/hasse_diagram.h"

namespace polymake { namespace fan{

   struct CellularDecoration : public GenericStruct<CellularDecoration> {
     DeclSTRUCT( DeclFIELD(face, Set<int>)
                 DeclFIELD(rank,int) );

     CellularDecoration() {}
     CellularDecoration(const Set<int>& f, int r) : face(f), rank(r) {}
   };
   
   class CustomDecorator {
      
      public:
         CustomDecorator(){}

         CellularDecoration compute_initial_decoration(const Set<int>& face) const {
            return CellularDecoration(face, 0);
         }

         CellularDecoration compute_decoration(const Set<int>& face, const CellularDecoration& bd) const {
            return CellularDecoration(face, bd.rank+1);
         }

         CellularDecoration compute_artificial_decoration(const NodeMap<Directed, CellularDecoration>& decor, const std::list<int>& max_faces) const {
            Set<int> D;
            D += -1;
            return CellularDecoration(D, -1);
         }
         
   };


} // namespace fan
} // namespace polymake
#endif