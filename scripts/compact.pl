use strict;
use warnings;

use application "fan";

my $k3 = load("k3.pcom");
$k3->COMPACTIFICATION;
# $k3 = new PolyhedralComplex(check_fan_objects(cube(3)));
# print $k3->COMPACTIFICATION->SEDENTARITY_SORT,"\n";
for(my $i=0; $i<3; $i++){
   my $f1 = $k3->compact_fcosheaf($i);
   $Polymake::User::Verbose::cpp = 3;
   my $d1 = build_full_chain($k3->COMPACTIFICATION, $k3->COMPACTIFICATION->ORIENTATIONS, $f1, false);
   # print $d1;
   print topaz::betti_numbers($d1),"\n";
}
# print $k3->COMPACTIFICATION->DECORATION;
