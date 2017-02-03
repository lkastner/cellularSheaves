application "fan";
; # W1-Sheaf on a 3-dim cube
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
print $pc->HASSE_DIAGRAM->FACES;
$w1 = $pc->wsheaf(1);
; #Building a sheaf from bases:
print $w1->BASES->type->full_name;
print $w1->BASES->{new Set<Int>([0,1,2,3,4,5,6,7])};
print $w1->BASES->{new Set<Int>([0,1,4,5])};
print $w1->BASES->{new Set<Int>([0,2])};
print $w1->BLOCKS->type->full_name;
print $w1->BLOCKS->{new Set<Set<Int>>([[0,1,2,3,4,5,6,7], [0,1,4,5]])};
print $w1->BLOCKS->{new Set<Set<Int>>([[0,1,4,5], [0,2]])};
; # Building the chain complex
print $pc->ORIENTATIONS->type->full_name;
print $pc->ORIENTATIONS->{new Set<Set<Int>>([[0,1,2,3,4,5,6,7], [0,1,4,5]])};
print $pc->ORIENTATIONS->{new Set<Set<Int>>([[0,1,4,5], [0,2]])};
$cs1 = $pc->compact_support_complex($w1);
$cs1->print();
print $cs1->BETTI_NUMBERS->type->full_name;
print $cs1->BETTI_NUMBERS;
print $cs1->COHOMOLOGIES->type->full_name;
print $cs1->COHOMOLOGIES;

; # Betti vs h-vector
@betti = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wsheaf($i);
   my $s = $pc->usual_cochain_complex($w);
   push @betti, $s->BETTI_NUMBERS;   
}
print new Matrix(@betti);

$cube = polytope::cube(3);
print $cube->DUAL_H_VECTOR;

; # Orlik-Solomon algebra of a matroid (tropical line)
$m = matroid::uniform_matroid(2,3);
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

; # Orlik-Solomon algebra of a matroid
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

; # Orlik-Solomon algebra of a valuated matroid
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

