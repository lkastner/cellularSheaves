use strict;
use warnings;

use application "topaz";
use application "fan";

my $k3 = load("k3.pcom");
$k3->COMPACTIFICATION;
# $k3 = new PolyhedralComplex(check_fan_objects(cube(3)));
# print $k3->COMPACTIFICATION->SEDENTARITY_SORT,"\n";
for(my $i=0; $i<3; $i++){
   my $f1 = $k3->compact_fcosheaf($i);
   verify_sheaf($f1->BLOCKS, $k3->COMPACTIFICATION->ADJACENCY,1);
   # save($f1, "sheaf$i");
   # $Polymake::User::Verbose::cpp = 3;
   my $d1 = new topaz::ChainComplex<Matrix<Rational>>(build_full_chain($k3->COMPACTIFICATION, $k3->COMPACTIFICATION->ORIENTATIONS, $f1, false));
   # # print $d1;
   my $betti = topaz::betti_numbers($d1);
   print $betti,"\n";
}
# print $k3->COMPACTIFICATION->DECORATION;
