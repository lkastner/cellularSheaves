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
print $t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t); 
$f0 = $berg->fcosheaf(0);
$f1 = $berg->fcosheaf(1);
$f2 = $berg->fcosheaf(2);
$s0 = $berg->usual_chain_complex($f0);
$s1 = $berg->usual_chain_complex($f1);
$s2 = $berg->usual_chain_complex($f2);
$s0->print();
$s1->print();
$s2->print();
print $s0->BETTI_NUMBERS;
print $s1->BETTI_NUMBERS;
print $s2->BETTI_NUMBERS;
$bm0 = $berg->borel_moore_complex($f0);
$bm1 = $berg->borel_moore_complex($f1);
$bm2 = $berg->borel_moore_complex($f2);
$bm0->print();
$bm1->print();
$bm2->print();
print $bm0->BETTI_NUMBERS;
print $bm1->BETTI_NUMBERS;
print $bm2->BETTI_NUMBERS;

application "graph";
$g = complete(4);
application "matroid";
$m = matroid_from_graph($g);
print $m->BASES;
$t = tropical::matroid_fan<Max>($m);
print $t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);
@result1 = ();
@result2 = ();
for(my $i=0; $i<3; $i++){
   my $f = $berg->fcosheaf($i);
   my $s = $berg->usual_chain_complex($f);
   my $bm = $berg->borel_moore_complex($f);
   push @result1, $s->BETTI_NUMBERS;
   push @result2, $bm->BETTI_NUMBERS;
}
print new Matrix(@result1);
print new Matrix(@result2);



$m = matroid::uniform_matroid(3,6);
$t = tropical::matroid_fan<Max>($m);
print $t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);

@result1 = ();
@result2 = ();
for(my $i=0; $i<6; $i++){
   my $f = $berg->fcosheaf($i);
   my $s = $berg->usual_chain_complex($f);
   my $bm = $berg->borel_moore_complex($f);
   push @result1, $s->BETTI_NUMBERS;
   push @result2, $bm->BETTI_NUMBERS;
}
print new Matrix(@result1);
print new Matrix(@result2);



## Tropical linear spaces 


$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
@result1 = ();
@result2 = ();
for(my $i=0;$i<3;$i++){
   my $fi = $tls->fcosheaf($i);
   my $si=$tls->usual_chain_complex($fi);
   my $bmi=$tls->borel_moore_complex($fi);
   push @result1, $si->BETTI_NUMBERS;
   push @result2, $bmi->BETTI_NUMBERS;
}  
print new Matrix(@result1);
print new Matrix(@result2);



### Some tropical hypersurfaces

application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = tropical::divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
print $div->F_VECTOR;
print $div->BOUNDED_FACES;

@result1 = ();
@result2 = ();
for(my $i=0;$i<3;$i++){
   my $fi = $div->fcosheaf($i);
   my $si=$div->usual_chain_complex($fi);
   my $bmi=$div->borel_moore_complex($fi);
   push @result1, $si->BETTI_NUMBERS;
   push @result2, $bmi->BETTI_NUMBERS;
} 
print new Matrix(@result1);
print new Matrix(@result2);



application "tropical";
$g = toTropicalPolynomial("max(0,x,y, 2*x - 2, 2*y-2, x+y-1, 3*x-6, 3*y-6, 2*x+y - 4,2*y + x- 4, 4*x -12, 4*y-12,  3*x+y -9,3*y + x- 9, 2*x + 2*y - 8)");

$gen3 = tropical::divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($g));  
application "fan";
$f0 = $gen3->fcosheaf(0);
$f1 = $gen3->fcosheaf(1);

$us0 = $gen3->usual_chain_complex($f0);
$us1 = $gen3->usual_chain_complex($f1);

$bm0 = $gen3->borel_moore_complex($f0);
$bm1 = $gen3->borel_moore_complex($f1);

$us0->print();
$us1->print();

print $us0->BETTI_NUMBERS;
print $us1->BETTI_NUMBERS;

$bm0->print();
$bm1->print();

print $bm0->BETTI_NUMBERS;
print $bm1->BETTI_NUMBERS;



application "tropical";
$f = toTropicalPolynomial("max(0,x,y,z, 2*x - 2, 2*y-2, 2*z-2, x+y-1, x+z-1, y+z-1, 3*x-6, 3*y-6, 3*z-6, 2*x+y - 4,2*y + x- 4, 2*x + z - 4, 2*z+x - 4, 2*y + z- 4, 2*z+y - 4, x+y+z +1, 4*x -12, 4*y-12, 4*z-12,  3*x+y -9,3*y + x- 9, 3*x + z - 9, 3*z+x - 9, 3*y + z- 9, 3*z+y - 9, 2*x + 2*y - 8, 2*x + 2*z - 8, 2*y + 2*z - 8,2*x + y +z - 7, x + 2*z +y - 7, 2*y + z + x- 7  )");

$k3 = tropical::divisor( (projective_torus<Max>(3)) , rational_fct_from_affine_numerator($f));

application "fan";
print $k3->BOUNDED_FACES->{0}->size;

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

