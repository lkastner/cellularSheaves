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

#include "polymake/fan/ChainComplexBuilder.h"

namespace polymake { namespace fan{
 
   typedef Lattice<BasicDecoration, lattice::Nonsequential> HasseDiagramType;

   class TrivialSelector {
      public:
         TrivialSelector(){}

         bool isValid(const Set<int>& face) const{
            return true;
         }
   };
   
   class NonFarSelector {
      private:
         Set<int> farFace;
      public:
         NonFarSelector(const Set<int>& ff) : farFace(ff) {}

         bool isValid(const Set<int>& face) const{
            return !(face - farFace).empty();
         }
   };
   
   class BoundedSelector {
      private:
         Set<int> farFace;
      public:
         BoundedSelector(const Set<int>& ff) : farFace(ff) {}

         bool isValid(const Set<int>& face) const{
            return (face * farFace).empty();
         }
   };

   topaz::ChainComplex<Matrix<Rational>> build_nonfar_chain(const HasseDiagramType& hd, const EdgeMap<Directed, int>& orientations, const Set<int>& ff, perl::Object cosheaf, bool cochain){
      NonFarSelector bs(ff);
      topaz::ChainComplex<Matrix<Rational>> result(build_chain_complex_from_hasse(hd, orientations, cosheaf, bs, cochain));
      return result;
   }

   
   topaz::ChainComplex<Matrix<Rational>> build_bounded_chain(const HasseDiagramType& hd, const EdgeMap<Directed, int>& orientations, const Set<int>& ff, perl::Object cosheaf, bool cochain){
      BoundedSelector bs(ff);
      topaz::ChainComplex<Matrix<Rational>> result(build_chain_complex_from_hasse(hd, orientations, cosheaf, bs, cochain));
      return result;
   }
   
   template<typename Decoration, typename SeqType>
   topaz::ChainComplex<Matrix<Rational>> build_full_chain(const Lattice<Decoration, SeqType>& hd, const EdgeMap<Directed, int>& orientations, perl::Object cosheaf, bool cochain){
      TrivialSelector ts;
      return build_chain_complex_from_hasse(hd, orientations, cosheaf, ts, cochain);
   }
   
   Function4perl(&build_bounded_chain, "build_bounded_chain( $ , $ , $ , $ , $ )");
   
   FunctionTemplate4perl("build_full_chain<Decoration, SeqType>( Lattice<Decoration, SeqType> , $ , $ , $ )");

   Function4perl(&build_nonfar_chain, "build_nonfar_chain( $ , $ , $ , $ , $ )");

} // namespace fan
} // namespace polymake
