## For theses examples h^{p,q} = 0 if p!=q. For h^{p,p} we recover 
## the h-vector of the polytope we took the cone over.

application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
print $pc->HASSE_DIAGRAM->FACES;
print $pc->ORIENTATIONS;
$w1 = $pc->wsheaf(1);
print $w1->BASES;
print $w1->BLOCKS;
print $w1->CHAIN_COMPLEX->BETTI_NUMBERS;
print $w1->CHAIN_COMPLEX->HOMOLOGIES;

$w2 = $pc->wsheaf(2);
print $w2->BASES;
print $w2->BLOCKS;
print $w2->CHAIN_COMPLEX->BETTI_NUMBERS;

$w3 = $pc->wsheaf(3);
print $w3->CHAIN_COMPLEX->BETTI_NUMBERS;

$w0 = $pc->wsheaf(0);
print $w0->CHAIN_COMPLEX->BETTI_NUMBERS;

@result = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wsheaf($i);
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
print new Matrix(@result);

$cube = polytope::cube(3);
print $cube->DUAL_H_VECTOR;

$pc = new PolyhedralComplex(check_fan_objects(new Cone(simplex(5))));
@result = ();
for(my $i=0; $i<6; $i++){
   my $w = $pc->wsheaf($i);
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
print new Matrix(@result);

$simplex = polytope::simplex(5);
print $simplex->DUAL_H_VECTOR;


# Bergman fans and tropical linear spaces


$m = matroid::uniform_matroid(2,3);
application "tropical";
$t = matroid_fan<Max>($m);
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
print $m->BASES;
application "tropical";
$t = matroid_fan<Max>($m);
application "fan";
$berg = new PolyhedralComplex($t);
$f0 = $berg->fcosheaf(0);
$f1 = $berg->fcosheaf(1);
$f2 = $berg->fcosheaf(2);
$f3 = $berg->fcosheaf(3);
$s0 = $berg->usual_chain_complex($f0);
$s1 = $berg->usual_chain_complex($f1);
$s2 = $berg->usual_chain_complex($f2);
$s3 = $berg->usual_chain_complex($f3);
print $s0->BETTI_NUMBERS;
print $s1->BETTI_NUMBERS;
print $s2->BETTI_NUMBERS;
print $s3->BETTI_NUMBERS;

$bm0 = $berg->borel_moore_complex($f0);
$bm1 = $berg->borel_moore_complex($f1);
$bm2 = $berg->borel_moore_complex($f2);
$bm3 = $berg->borel_moore_complex($f3);
print $bm0->BETTI_NUMBERS;
print $bm1->BETTI_NUMBERS;
print $bm2->BETTI_NUMBERS;
print $bm3->BETTI_NUMBERS;


$m = matroid::uniform_matroid(3,6);
$t = matroid_fan<Max>($m);
application "fan";
$berg = new PolyhedralComplex($t);

for(my $i=0; $i<6; $i++){
   my $f = $berg->fcosheaf($i);
   my $bm = $berg->usual_chain_complex($f);
   push @result, $bm->BETTI_NUMBERS;
}
print new Matrix(@result);


@result = ();
for(my $i=0; $i<6; $i++){
   my $f = $berg->fcosheaf($i);
   my $bm = $berg->borel_moore_complex($f);
   push @result, $bm->BETTI_NUMBERS;
}
print new Matrix(@result);


## Tropical linear spaces 


$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
@result = ();
for(my $i=0;$i<3;$i++){
   my $fi = $tls->fsheaf($i);
   my $si=$tls->usual_chain_complex($fi);
   push @result, $si->BETTI_NUMBERS;
}  
print new Matrix(@result);

@result = ();
for(my $i=0;$i<3;$i++){
   my $fi = $tls->fsheaf($i);
   my $si=$tls->borel_moore_complex($fi);
   push @result, $si->BETTI_NUMBERS;
}  
print new Matrix(@result);


### Some tropical hypersurfaces

application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = tropical::divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
print $div->F_VECTOR;
print $div->BOUNDED_FACES;

@result = ();
for(my $i=0;$i<3;$i++){
   my $fi = $div->fsheaf($i);
   my $si=$div->usual_chain_complex($fi);
   push @result, $si->BETTI_NUMBERS;
} 
print new Matrix(@result);

@result = ();
for(my $i=0;$i<3;$i++){
   my $fi = $div->fsheaf($i);
   my $si=$div->borel_moore_complex($fi);
   push @result, $si->BETTI_NUMBERS;
} 
print new Matrix(@result);







 

