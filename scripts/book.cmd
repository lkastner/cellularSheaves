application "fan";
$pc = new PolyhedralComplex( 
   check_fan_objects(new Cone(cube(3))));
@betti = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wsheaf($i);
   my $s = $pc->usual_cochain_complex($w);
   push @betti, $s->BETTI_NUMBERS;   
}
print new Matrix(@betti);

$cube = polytope::cube(3);
print $cube->DUAL_H_VECTOR;


application "matroid";
$m = uniform_matroid(2,3);
application "tropical";
$t = matroid_fan<Max>($m);
$t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);
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

application "graph";
$g = complete(4);
application "matroid";
$m = matroid_from_graph($g);
application "tropical";
$t = matroid_fan<Max>($m);
$t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);
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
$t = tropical::matroid_fan<Max>($m);
$t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);
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

@wbetti_usual = ();
@wbetti_cs = ();
for(my $i=0;$i<3;$i++){
   my $wi = $tls->wsheaf($i);
   my $wsi=$tls->usual_cochain_complex($wi);
   my $wcsi=$tls->compact_support_complex($wi);
   push @wbetti_usual, $wsi->BETTI_NUMBERS;
   push @wbetti_cs, $wcsi->BETTI_NUMBERS;
}
print new Matrix(@wbetti_usual);
print new Matrix(@wbetti_cs);


application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = divisor( (projective_torus<Max>(2) ),
   rational_fct_from_affine_numerator($f));
application "fan";
@betti_usual = ();
@betti_bm = ();
for(my $i=0;$i<2;$i++){
   my $fi = $div->fcosheaf($i);
   my $si=$div->usual_chain_complex($fi);
   my $bmi=$div->borel_moore_complex($fi);
   push @betti_usual, $si->BETTI_NUMBERS;
   push @betti_bm, $bmi->BETTI_NUMBERS;
}
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);


application "tropical";
$f = toTropicalPolynomial("max(0,x,y,z, 2*x-2, 
   2*y-2, 2*z-2, x+y-1, x+z-1, y+z-1, 3*x-6, 
   3*y-6, 3*z-6, 2*x+y-4, 2*y+x-4, 2*x+z-4, 
   2*z+x-4, 2*y+z-4, 2*z+y-4, x+y+z+1, 4*x-12, 
   4*y-12, 4*z-12, 3*x+y-9, 3*y+x-9, 3*x+z-9, 
   3*z+x-9, 3*y+z-9, 3*z+y-9, 2*x+2*y-8, 
   2*x+2*z-8, 2*y+2*z-8, 2*x+y+z-7, x+2*z+y-7, 
   2*y+z+x-7)");
$k3 = divisor((projective_torus<Max>(3)),
   rational_fct_from_affine_numerator($f));
application "fan";
@betti_usual = ();
@betti_bm = ();
for(my $i=0;$i<3;$i++){
   my $fi = $k3->fcosheaf($i);
   my $si=$k3->usual_chain_complex($fi);
   my $bmi=$k3->borel_moore_complex($fi);
   push @betti_usual, $si->BETTI_NUMBERS;
   push @betti_bm, $bmi->BETTI_NUMBERS;
}
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);


application "tropical";
$f = toTropicalPolynomial("max(0,x,y,z, 2*x-2, 
   2*y-2, 2*z-2, x+y-1, x+z-1, y+z-1, 3*x-6, 
   3*y-6, 3*z-6, 2*x+y-4, 2*y+x-4, 2*x+z-4, 
   2*z+x-4, 2*y+z-4, 2*z+y-4, x+y+z+1, 4*x-12, 
   4*y-12, 4*z-12, 3*x+y-9, 3*y+x-9, 3*x+z-9, 
   3*z+x-9, 3*y+z-9, 3*z+y-9, 2*x+2*y-8, 
   2*x+2*z-8, 2*y+2*z-8, 2*x+y+z-7, x+2*z+y-7, 
   2*y+z+x-7)");
$k3 = divisor((projective_torus<Max>(3)),
   rational_fct_from_affine_numerator($f));
application "fan";
@numbers = (0..2);
@cosheaves = map{$k3->fcosheaf($_)} @numbers;
@usualChainComplexes = map{$k3->usual_chain_complex($_)} @cosheaves;
@bmComplexes = map{$k3->borel_moore_complex($_)} @cosheaves;
@betti_usual = map{$_->BETTI_NUMBERS} @usualChainComplexes;
@betti_bm = map{$_->BETTI_NUMBERS} @bmComplexes;
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);
$usualChainComplexes[0]->print();


$f0 = $k3->fcosheaf(0);
$f1 = $k3->fcosheaf(1);
$f2 = $k3->fcosheaf(2);

$us0 = $k3->usual_chain_complex($f0);
$us1 = $k3->usual_chain_complex($f1);
$us2 = $k3->usual_chain_complex($f2);

$bm0 = $k3->borel_moore_complex($f0);
$bm1 = $k3->borel_moore_complex($f1);
$bm2 = $k3->borel_moore_complex($f2);


$us0->print();
print $us0->BETTI_NUMBERS;
$us1->print();
print $us1->BETTI_NUMBERS;
$us2->print();
print $us2->BETTI_NUMBERS;
$bm0->print();
print $bm0->BETTI_NUMBERS;
$bm1->print();
print $bm1->BETTI_NUMBERS;
$bm2->print();
print $bm2->BETTI_NUMBERS;

