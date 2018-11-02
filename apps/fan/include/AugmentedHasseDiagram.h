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

#ifndef AUGMENTED_HASSE_DIAGRAM
#define AUGMENTED_HASSE_DIAGRAM

#include "polymake/client.h"
#include "polymake/GenericMatrix.h"
#include "polymake/Map.h"
#include "polymake/Set.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"
#include "polymake/fan/CellularData.h"
#include "polymake/fan/CellularClosure.h"
#include "polymake/fan/CellularDecoration.h"

namespace polymake { namespace fan{
  
   template<typename DecoratorType>
   class AugmentedHasseDiagram {

      private:
         const CellularData& cd;
         const DecoratorType& decorator;
         using DecorationType=typename DecoratorType::DecorationType;
         CellularClosureOperator<DecorationType> tco;
         Lattice<DecorationType, lattice::Nonsequential> hasseDiagram;

         void compute_hasse_diagram() {
            hasseDiagram = graph::lattice_builder::compute_lattice_from_closure<DecorationType>(tco, TrivialCut<DecorationType>(), decorator, true, std::false_type());
         }

      public:
         AugmentedHasseDiagram(const CellularData& cd_in, DecoratorType& decorator_in, perl::Object pc) : cd(cd_in), decorator(decorator_in), tco(cd, pc) {
            compute_hasse_diagram();
         }

         void print() const {
            cout << "Vertices of compactification:" << endl;
            for(const auto& v:cd.int2vertices){
               cout << v.first << ": " << v.second << endl;
            }
            cout << endl;
            cout << "Hasse diagram of compactification:" << endl;
            const auto& G = hasseDiagram.graph();
            const auto& D = hasseDiagram.decoration();
            for(const auto& node: nodes(G)){
               cout << node << ": ";
               cout << "neighbors: " << G.out_adjacent_nodes(node) << " ";
               cout << "decoration: " << D[node] << endl;
            }
         }

         const CellularClosureOperator<DecorationType>& get_ClosureOperator(){
            return tco;
         }

         const Lattice<DecorationType, lattice::Nonsequential>& get_HasseDiagram() const {
            return hasseDiagram;
         }

   };

} // namespace fan
} // namespace polymake

#endif
