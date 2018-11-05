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

namespace polymake { namespace fan{
  
 



   // perl::Object
   void tropcomp(perl::Object pc){
      CompactificationData cd(pc);
      CellularDecorator decorator(cd);
      AugmentedHasseDiagram<CellularDecorator> AHD(cd, decorator, pc);
      AHD.print();

      // Use a different decorator:
      cout << endl << "-------------------------------" << endl;
      cout << "Using SedentarityDecorator: " << endl;
      SedentarityDecorator sd(cd);
      AugmentedHasseDiagram<SedentarityDecorator> AHDSD(cd, sd, pc);
      AHDSD.print();
    
      // const Lattice<CellularDecoration, lattice::Nonsequential>& HD(AHD.get_HasseDiagram());
      // cout << "Graph: " << endl;
      // cout << HD.graph() << endl;
      // cout << "Node decoration: " << endl;
      // cout << HD.decoration() << endl;
      // return HD.makeObject();
      
     //  template <typename Decoration, typename ClosureOperator, typename CrossCut, typename Decorator, bool dual, typename SeqType = lattice::Nonsequential >
     //  Lattice<Decoration, SeqType> compute_lattice_from_closure(
     //        ClosureOperator cl,
     //        const CrossCut& cut,
     //        const Decorator& decorator,
     //        bool wants_artificial_top_node,
     //        bool_constant<dual> built_dually,
     //        Lattice<Decoration, SeqType> lattice = Lattice<Decoration>(),
     //        Set<int> queuing_nodes = Set<int>())
   }

   Function4perl(&tropcomp, "tropcomp( $ )");

} // namespace fan
} // namespace polymake
