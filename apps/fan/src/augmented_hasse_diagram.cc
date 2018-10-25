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
#include "polymake/fan/hasse_diagram.h"

namespace polymake { namespace fan{

/*
application "fan";
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[1,1,0],[1,0,1],[0,1,1]], INPUT_CONES=>[[0,1,3],[1,2,3]]);


*/
  
	using graph::Lattice;
   using namespace graph::lattice;
   using namespace fan::lattice;


   struct IteratorWrap {
      FacetList groundSet;
      FacetList::const_iterator state;

      IteratorWrap(FacetList&& gs): groundSet(gs){
         state = groundSet.begin();
      }

      Set<int> operator*(){
         return *state;
      }

      void operator++(){
         state++;
      }

      bool at_end(){
         return state == groundSet.end();
      }

   };

 
   class AugmentedHasseDiagram {

      private:
           FaceMap<> face_index_map;

         Map<int, Set<int>> int2vertices;
         Map<Set<int>, int> vertices2int;
         int nVertices;
         Set<int> farVertices;
         ComplexPrimalClosure<> oldClosureOperator;

         IncidenceMatrix<> construct_old_closure_operator(perl::Object pc) {
            RestrictedIncidenceMatrix<> building_matrix;
            const Array<IncidenceMatrix<> >& maximal_vifs = pc.give("MAXIMAL_CONES_INCIDENCES");
            bool is_pure, is_complete = false;
            pc.give("PURE") >> is_pure;
            TopologicalType tt(is_pure, is_complete);
            const IncidenceMatrix<>& maximal_cones = pc.give("MAXIMAL_CONES");
            const int n_vertices = maximal_cones.cols();
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
         typedef Set<int> ClosureData;
         AugmentedHasseDiagram(perl::Object pc) : oldClosureOperator(construct_old_closure_operator(pc)) {
            pc.give("FAR_VERTICES") >> farVertices;
            const Matrix<Rational>& vertices = pc.give("VERTICES");
            nVertices = vertices.rows();
            const Lattice<BasicDecoration, Nonsequential>& oldHasseDiagram = pc.give("HASSE_DIAGRAM");
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

         Set<int> old_closure(const Set<int>& a) const {
            Set<int> result(range(0,nVertices));
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

         Set<int> closure(const Set<int> a) const {
            Set<int> originalRealisation;
            for(const auto i:a){
               originalRealisation += int2vertices[i];
            }
            Set<int> originalClosure = old_closure(originalRealisation);
            Set<int> commonRays = originalRealisation * farVertices;
            for(const auto i : a){
               commonRays = commonRays * int2vertices[i];
            }
            Set<int> result;
            for(const auto& v:vertices2int){
               if(incl(commonRays, v.first)<=0 && incl(v.first, originalClosure)<=0){
                  result += v.second;
               }
            }
            return result;
         }

         void print() const {
            for(const auto& v:int2vertices){
               cout << v.first << ": " << v.second << endl;
            }
         }
         
         Set<int> closure_of_empty_set(){
            Set<int> empty;
            return empty;
         }

         FaceIndexingData get_indexing_data(const ClosureData& data)
        {
          int& fi = face_index_map[data];
          return FaceIndexingData(fi, fi == -1, fi == -2);
        }

        Set<int> compute_closure_data(const BasicDecoration& bd) const {
            return bd.face;
        }

        IteratorWrap get_closure_iterator(const Set<int>& face) const {
          Set<int> all = pm::range(0,int2vertices.size()-1);
          Set<int> toadd = all-face;
          FacetList result;
          for(auto i:toadd){
            result.insertMin(closure(face+i));
         }
         return IteratorWrap(std::move(result));
        }


   };


   class TrivialDecorator {
      
      public:
         TrivialDecorator(){}

         BasicDecoration compute_initial_decoration(const Set<int>& face) const {
            return BasicDecoration(face, 0);
         }

         BasicDecoration compute_decoration(const Set<int>& face, const BasicDecoration& bd) const {
            return BasicDecoration(face, bd.rank+1);
         }

         BasicDecoration compute_artificial_decoration(const NodeMap<Directed, BasicDecoration>& decor, const std::list<int>& max_faces) const {
            Set<int> D;
            D += -1;
            return BasicDecoration(D, -1);
         }
         
   };

   perl::Object tropcomp(perl::Object pc){
      AugmentedHasseDiagram AHD(pc);
      AHD.print();
      Set<int> S, S1;
      S += 8;
      S += 7;
      S1 += 7;
      cout << "S: " << S << " S1: " << S1 << endl;
      cout << "incl(S1, S) " << incl(S1, S) << endl;
      cout << "incl(S, S) " << incl(S, S) << endl;
      cout << "incl(S, S1) " << incl(S, S1) << endl;
      cout << AHD.closure(S) << endl;
      S.clear();
      S += 1;
      S += 2;
      cout << AHD.closure(S) << endl;
    
      TrivialDecorator decorator;
      Lattice<BasicDecoration, lattice::Nonsequential> HD=graph::lattice_builder::compute_lattice_from_closure<BasicDecoration>(AHD, TrivialCut<BasicDecoration>(), decorator, true, std::false_type()); 
      cout << HD.graph() << endl;
      cout << HD.decoration() << endl;
      return HD.makeObject();
      
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
#endif
