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

#ifndef COMPACTIFICATION_DATA
#define COMPACTIFICATION_DATA

#include "polymake/client.h"
#include "polymake/graph/Decoration.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/BasicLatticeTypes.h"
#include "polymake/graph/lattice_builder.h"
#include "polymake/graph/LatticePermutation.h"
#include "polymake/fan/hasse_diagram.h"
#include "polymake/GenericMatrix.h"
#include "polymake/Map.h"
#include "polymake/Set.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"
namespace polymake { namespace fan{
   
   using graph::Lattice;
   using namespace graph::lattice;
   using namespace fan::lattice;

   class CompactificationData {
      public:
         Map<Int, Set<Int>> int2vertices;
         Map<Set<Int>, Int> vertices2int;
         Int nVertices;
         Set<Int> farVertices;
         Matrix<Rational> vertices;

      public:
         CompactificationData(perl::Object pc){
            pc.give("FAR_VERTICES") >> farVertices;
            pc.give("VERTICES") >> vertices;
            nVertices = vertices.rows();
            Set<Int> topNode; topNode += -1;
            // cout << oldHasseDiagram << endl;
            Int i = 0;
            // Build new vertices
            const Lattice<BasicDecoration, Nonsequential>& oldHasseDiagram(pc.give("HASSE_DIAGRAM"));
            for(const auto& f : oldHasseDiagram.decoration()){
               if(f.face != topNode) { 
                  Int faceDim = f.rank-1;
                  Int tailDim = rank(vertices.minor(f.face * farVertices, All));
                  if(faceDim == tailDim){
                     int2vertices[i] = f.face;
                     vertices2int[f.face] = i;
                     i++;
                  } else {
                  }
               }
            }
         }

   };

} // namespace fan
} // namespace polymake

#endif
