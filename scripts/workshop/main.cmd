application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
print $pc->HASSE_DIAGRAM->FACES;
print $pc->ORIENTATIONS->{new Set<Set<Int>>([[0,1,2,3,4,5,6,7], [0,1,4,5]])};
print $pc->ORIENTATIONS->{new Set<Set<Int>>([[0,1,4,5], [0,2]])};
$w1 = $pc->wsheaf(1);
print $w1->BASES->{new Set<Int>([0,1,2,3,4,5,6,7])};
print $w1->BASES->{new Set<Int>([0,1,4,5])};
print $w1->BASES->{new Set<Int>([0,2])};
print $w1->BLOCKS->{new Set<Set<Int>>([[0,1,2,3,4,5,6,7], [0,1,4,5]])};
print $w1->BLOCKS->{new Set<Set<Int>>([[0,1,4,5], [0,2]])};
$cs1 = $pc->compact_support_complex($w1);
$cs1->print();
print $cs1->BETTI_NUMBERS;
print $cs1->COHOMOLOGIES;

@betti = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wsheaf($i);
   my $s = $pc->usual_cochain_complex($w);
   push @betti, $s->BETTI_NUMBERS;   
}
print new Matrix(@betti);

$cube = polytope::cube(3);
print $cube->DUAL_H_VECTOR;


$m = matroid::uniform_matroid(2,3);
$berg = tropical::matroid_fan<Max>($m);
$f0 = $berg->fcosheaf(0);
$f1 = $berg->fcosheaf(1);
$f2 = $berg->fcosheaf(2);
$s0 = $berg->usual_chain_complex($f0);
$s1 = $berg->usual_chain_complex($f1);
$s2 = $berg->usual_chain_complex($f2);
print $s0->BETTI_NUMBERS;
print $s1->BETTI_NUMBERS;
print $s2->BETTI_NUMBERS;

$bm0 = $berg->borel_moore_complex($f0);
$bm1 = $berg->borel_moore_complex($f1);
$bm2 = $berg->borel_moore_complex($f2);
print $bm0->BETTI_NUMBERS;
print $bm1->BETTI_NUMBERS;
print $bm2->BETTI_NUMBERS;

$g = graph::complete(4);
$m = matroid::matroid_from_graph($g);
$berg = tropical::matroid_fan<Max>($m);
@betti_usual = ();
@betti_bm = ();
for(my $i=0; $i<3; $i++){
   my $f = $berg->fcosheaf($i);
   my $s = $berg->usual_chain_complex($f);
   my $bm = $berg->borel_moore_complex($f);
   push @betti_usual, $s->BETTI_NUMBERS;
   push @betti_bm, $bm->BETTI_NUMBERS;
}
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);


$m = matroid::uniform_matroid(3,6);
$berg = tropical::matroid_fan<Max>($m);
@betti_usual = ();
@betti_bm = ();
for(my $i=0; $i<3; $i++){
   my $f = $berg->fcosheaf($i);
   my $s = $berg->usual_chain_complex($f);
   my $bm = $berg->borel_moore_complex($f);
   push @betti_usual, $s->BETTI_NUMBERS;
   push @betti_bm, $bm->BETTI_NUMBERS;
}
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);


$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>( 
   BASES=>matroid::uniform_matroid(3,6)->BASES, 
   VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
@betti_usual = ();
@betti_bm = ();
for(my $i=0;$i<3;$i++){
   my $fi = $tls->fcosheaf($i);
   my $si=$tls->usual_chain_complex($fi);
   my $bmi=$tls->borel_moore_complex($fi);
   push @betti_usual, $si->BETTI_NUMBERS;
   push @betti_bm, $bmi->BETTI_NUMBERS;
}
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);

load_commands("wsheaf.cmd");
load_commands("tropical.cmd");

