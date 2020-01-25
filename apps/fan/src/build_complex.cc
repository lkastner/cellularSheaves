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

#include "polymake/fan/ChainComplexBuilder.h"

namespace polymake { namespace fan{
 
   typedef Lattice<BasicDecoration, lattice::Nonsequential> HasseDiagramType;

   class TrivialSelector {
      public:
         TrivialSelector(){}

         bool isValid(const Set<Int>& face) const{
            return true;
         }
   };
   
   class NonFarSelector {
      private:
         Set<Int> farFace;
      public:
         NonFarSelector(const Set<Int>& ff) : farFace(ff) {}

         bool isValid(const Set<Int>& face) const{
            return !(face - farFace).empty();
         }
   };
   
   class BoundedSelector {
      private:
         Set<Int> farFace;
      public:
         BoundedSelector(const Set<Int>& ff) : farFace(ff) {}

         bool isValid(const Set<Int>& face) const{
            return (face * farFace).empty();
         }
   };

   topaz::ChainComplex<Matrix<Rational>> build_nonfar_chain(const HasseDiagramType& hd, const EdgeMap<Directed, Int>& orientations, const Set<Int>& ff, BigObject cosheaf, bool cochain){
      NonFarSelector bs(ff);
      EdgeMap<Directed, Matrix<Rational>> blocks;
      cosheaf.give("BLOCKS") >> blocks;
      topaz::ChainComplex<Matrix<Rational>> result(build_chain_complex_from_hasse(hd, orientations, blocks, bs, cochain));
      return result;
   }

   
   topaz::ChainComplex<Matrix<Rational>> build_bounded_chain(const HasseDiagramType& hd, const EdgeMap<Directed, Int>& orientations, const Set<Int>& ff, BigObject cosheaf, bool cochain){
      BoundedSelector bs(ff);
      EdgeMap<Directed, Matrix<Rational>> blocks;
      cosheaf.give("BLOCKS") >> blocks;
      topaz::ChainComplex<Matrix<Rational>> result(build_chain_complex_from_hasse(hd, orientations, blocks, bs, cochain));
      return result;
   }
   
   template<typename Decoration, typename SeqType, typename Scalar>
   topaz::ChainComplex<Matrix<Scalar>> build_full_chain(const Lattice<Decoration, SeqType>& hd, const EdgeMap<Directed, Int>& orientations, const EdgeMap<Directed, Matrix<Scalar>>& blocks, bool cochain){
      TrivialSelector ts;
      return build_chain_complex_from_hasse(hd, orientations, blocks, ts, cochain);
   }
   
   Function4perl(&build_bounded_chain, "build_bounded_chain( $ , $ , $ , $ , $ )");
   
   FunctionTemplate4perl("build_full_chain<Decoration, SeqType, Scalar>( Lattice<Decoration, SeqType> , EdgeMap<Directed,Int> , EdgeMap<Directed, Matrix<Scalar>> , $ )");

   Function4perl(&build_nonfar_chain, "build_nonfar_chain( $ , $ , $ , $ , $ )");

} // namespace fan
} // namespace polymake
