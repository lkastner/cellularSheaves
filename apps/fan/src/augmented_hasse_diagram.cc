/* Copyright (c) 1997-2014
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

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

#ifndef POLYMAKE_WEDGE_MATRIX
#define POLYMAKE_WEDGE_MATRIX

#include "polymake/client.h"
#include "polymake/GenericMatrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/Map.h"
#include "polymake/Set.h"
#include "polymake/Integer.h"
#include "polymake/linalg.h"
#include "polymake/Array.h"
#include "polymake/PowerSet.h"
#include <list>
#include "polymake/graph/Decoration.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/BasicLatticeTypes.h"
#include "polymake/graph/lattice_builder.h"
#include "polymake/graph/LatticePermutation.h"

namespace polymake { namespace fan{
  
	using graph::Lattice;
   using graph::lattice::BasicDecoration;
   using graph::lattice::BasicDecorator;
   using graph::lattice::BasicClosureOperator;
   using graph::lattice::TrivialCut;
   using graph::lattice::SetAvoidingCut;
   using graph::lattice::RankCut;
   using graph::lattice::CutAnd;
   using graph::lattice::Sequential;
   using graph::lattice::Nonsequential;

 
   class AugmentedHasseDiagram {
      public:
         Matrix<Rational> vertices;
         Matrix<Rational> rays;
			

   };


} // namespace fan
} // namespace polymake
#endif
