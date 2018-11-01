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
#include "polymake/fan/CellularClosure.h"
#include "polymake/fan/CellularDecoration.h"

namespace polymake { namespace fan{
   
 
   class AugmentedHasseDiagram {

      private:
         Map<int, Set<int>> int2vertices;
         Map<Set<int>, int> vertices2int;
         int nVertices;
         Set<int> farVertices;
         CellularClosureOperator tco;
         Lattice<CellularDecoration, lattice::Nonsequential> hasseDiagram;

      public:
         AugmentedHasseDiagram(perl::Object pc) : tco(int2vertices, vertices2int, nVertices, farVertices, pc) {
            pc.give("FAR_VERTICES") >> farVertices;
            const Matrix<Rational>& vertices = pc.give("VERTICES");
            nVertices = vertices.rows();
            const Lattice<BasicDecoration, Nonsequential>& oldHasseDiagram(pc.give("HASSE_DIAGRAM"));
            Set<int> topNode; topNode += -1;
            // cout << oldHasseDiagram << endl;
            int i = 0;
            for(const auto& f : oldHasseDiagram.decoration()){
               if(f.face != topNode) { 
                  int faceDim = f.rank-1;
                  int tailDim = rank(vertices.minor(f.face * farVertices, All));
                  if(faceDim == tailDim){
                     int2vertices[i] = f.face;
                     vertices2int[f.face] = i;
                     i++;
                  } else {
                  }
               }
            }
         }

         void print() const {
            cout << "Vertices:" << endl;
            for(const auto& v:int2vertices){
               cout << v.first << ": " << v.second << endl;
            }
         }

         const CellularClosureOperator& get_ClosureOperator(){
            return tco;
         }

         void compute_hasse_diagram() {
            CustomDecorator decorator;
            hasseDiagram = graph::lattice_builder::compute_lattice_from_closure<CellularDecoration>(tco, TrivialCut<CellularDecoration>(), decorator, true, std::false_type());
         }

         const Lattice<CellularDecoration, lattice::Nonsequential>& get_HasseDiagram() const {
            return hasseDiagram;
         }

   };



   // perl::Object
   void tropcomp(perl::Object pc){
      AugmentedHasseDiagram AHD(pc);
      AHD.print();
      AHD.compute_hasse_diagram();
    
      const Lattice<CellularDecoration, lattice::Nonsequential>& HD(AHD.get_HasseDiagram());
      cout << "Graph: " << endl;
      cout << HD.graph() << endl;
      cout << "Node decoration: " << endl;
      cout << HD.decoration() << endl;
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
