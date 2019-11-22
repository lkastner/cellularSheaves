use strict;
use warnings;
use application "fan";

my $pc = new PolyhedralComplex(check_fan_objects(cube(3)));

# for(my $i = 0; $i<=3; $i++){
my $i = 1;

my $f1 = $pc->fcosheaf($i);
my $d1 = $pc->borel_moore_complex($f1);

my $cf1 = $pc->compact_fcosheaf($i);
my $comp = $pc->COMPACTIFICATION;
my $cdiff1 = build_full_chain($comp, $comp->ORIENTATIONS, $cf1, false);
my $cd1 = new ChainComplex(INPUT_DIFFERENTIALS=>$cdiff1);

# print $d1->BETTI_NUMBERS,"\n";
# print $cd1->BETTI_NUMBERS,"\n";
# 
# print $d1->INPUT_DIFFERENTIALS;
# print "-----\n";
# print $cdiff1;
# 
# print "-----\n";
for (my $e=entire(edges($pc->HASSE_DIAGRAM->ADJACENCY)); $e; ++$e) {
print $e->from_node," ",$e->to_node,": ",$f1->BLOCKS->[$$e];
}
print "-----\n";
print $cf1->BLOCKS;
# 
# print rows_labeled($pc->HASSE_DIAGRAM->DECORATION);
# print rows_labeled($comp->DECORATION);
# 
# print $cf1->BLOCKS->edge(0,1);
# print $f1->BLOCKS->edge(10,9);
print "---\n";
print $d1->BETTI_NUMBERS,"\n";
print $cd1->BETTI_NUMBERS,"\n";
# }
