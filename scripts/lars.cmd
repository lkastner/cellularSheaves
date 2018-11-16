application "fan";
$f = new PolyhedralFan(RAYS=>[[1,0,0],[1,1,0],[1,0,1],[1,1,1]], MAXIMAL_CONES=>[[0,1,2],[1,2,3]]);
$pc = new PolyhedralComplex($f);
print keys %{$pc->CHOSEN_BASES};

$s = new Set<Int>(1,2,3);
print $pc->CHOSEN_BASES->{$s};

$h = new Set<Int>([0]);
print $pc->CHOSEN_BASES->{$h}->type->full_name;


$t = new Set<Int>(2,3);
$pst = new Pair<Set<Int>, Set<Int> >($s, $t);
print $pc->SIMPLE_BLOCKS->{$pst};
print $pc->CHOSEN_BASES->{$s};
print $pc->CHOSEN_BASES->{$t};


#############################################
application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(4))));

$A = zero_matrix(4,3);
print choose_basis($A)->rows;
print choose_basis($A)->cols;


$A = new Matrix([[2, 1/2, 0],[1,0,0]]);
print dense(to_lattice_basis($A));


application "fan";
$A = new Matrix([[1,2,3,4],[3,4,5,6],[7,8,9,0]]);
print wedge_matrix($A, 2);

application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
@result = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wcosheaf($i);
   my $cs = $pc->borel_moore_complex($w);
   push @result, $cs->BETTI_NUMBERS;
}
print new Matrix(@result);

application "fan";
$d = 4;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$w2 = $pc->wsheaf(2);
$s = new Set<Int>(12,13,14,15);
$t = new Set<Int>(14,15);
$p = new Pair<Set<Int>, Set<Int> >($s, $t);
print $w2->BLOCKS->{$p};
print wedge_matrix($pc->SIMPLE_BLOCKS->{$p},2)->cols;


for(my $i = 1; $i<$d; $i++){
   my $wi = $pc->wcosheaf($i);
   my $cs = $pc->borel_moore_complex($wi);
   print $i,": ",$cs->BETTI_NUMBERS,"\n";
}

################################################################
# Testing c++ methods

application "fan";
$m = new Matrix([[1,2,3,4],[3,4,5,6],[5,6,7,8]]);
print wedge_matrix($m, 1);
print choose_basis($m);

application "fan";
$d = 4;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$s = $pc->HASSE_DIAGRAM->nodes_of_dim(3);
@s = map($pc->HASSE_DIAGRAM->FACES->[$_], @$s);
$sigmas = new Array<Set<Int> >(\@s);
$t = $pc->HASSE_DIAGRAM->nodes_of_dim(2);
@t = map($pc->HASSE_DIAGRAM->FACES->[$_], @$t);
$taus = new Array<Set<Int> >(\@t);
$orientations = $pc->ORIENTATIONS;
$blocks = $pc->SIMPLE_BLOCKS;
print assemble_matrix($sigmas, $taus, $blocks, $orientations);


###########################################################

#Testing f blocks

borel_moore_complex($f1);

print build_generator_matrix_f($tau, $pc, 1); 

###########################################################

application "graph";
$g = complete(4);
application "matroid";
$m = matroid_from_graph($g);
print $m->BASES;
application "tropical";
$t = matroid_fan<Max>($m);
application "fan";
print $t->RAYS;
$pcFan = new PolyhedralComplex($t);
$f1 = $pcFan->fcosheaf(1);
$boundedChain = build_chain_complex($f1->BLOCKS, $pcFan->BOUNDED_FACES, $pcFan->ORIENTATIONS);



###############################################################################


application "tropical"; 
$f = toTropicalPolynomial("max(0,x,y,z)"); 
$div = divisor( (projective_torus<Max>(3)) , rational_fct_from_affine_numerator($f)); 
application "fan"; 
$f0 = $div->fcosheaf(0); 
$f1 = $div->fcosheaf(1); 
$f2 = $div->fcosheaf(2);
$f3 = $div->fcosheaf(3); 
$bm0 = $div->borel_moore_complex($f0); 
$bm1 = $div->borel_moore_complex($f1); 
$bm2 = $div->borel_moore_complex($f2); 
$bm3 = $div->borel_moore_complex($f3); 
$bm0->print();
$bm1->print();
$bm2->print();
$bm3->print();
print $bm0->BETTI_NUMBERS;
print $bm1->BETTI_NUMBERS;
print $bm2->BETTI_NUMBERS;
print $bm3->BETTI_NUMBERS;


$us0 = $div->usual_chain_complex($f0);
$us1 = $div->usual_chain_complex($f1); 
$us2 = $div->usual_chain_complex($f2);
$us3 = $div->usual_chain_complex($f3); 

$us0->print();
$us1->print();
$us2->print();
$us3->print();

print $us0->BETTI_NUMBERS;
print $us1->BETTI_NUMBERS;
print $us2->BETTI_NUMBERS;
print $us3->BETTI_NUMBERS;

print $bm1->IS_WELLDEFINED;

print $div->ORIENTATIONS;
$u1 = $div->usual_chain_complex($f1);

$bm1 = $div->borel_moore_complex($f1); 
$bm1->print();
###############################################################################

application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = tropical::divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
print $div->INTERNAL_BASES;
print $div->ORIENTATIONS;
$f1 = $div->fcosheaf(1);


###############################################################################
application "fan"; 
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
print $pc->ORIENTATIONS->{new Set<Set<Int> >([[2,3],[3]])};
print $pc->ORIENTATIONS->{new Set<Set<Int> >([[2,3],[2]])};
$w0 = $pc->wsheaf(0);
$w1 = $pc->wsheaf(1);
$w2 = $pc->wsheaf(2);
$w3 = $pc->wsheaf(3);

###############################################################################
application "fan";
$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
$f1 = $tls->fcosheaf(1);
$s1=$tls->usual_chain_complex($f1);
$bm1=$tls->borel_moore_complex($f1);


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

###############################################################################

application "fan";
$d = 4;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$w2 = $pc->wsheaf(3);
$c = $pc->usual_cochain_complex($w2);
$cc = $c->INTERNAL_COMPLEX->INNER;
print topaz::homology<Integer>($cc,1);


#######################
# Compactification
application "fan";
$pc = new PolyhedralComplex(POINTS=>[[1,0,0],[0,1,0],[0,0,1]],INPUT_POLYTOPES=>[[0,1,2]]);
print $pc->COMPACTIFICATION->ALL_FACES;
print $pc->COMPACTIFICATION->VERTICES;
$O = $pc->COMPACTIFICATION->VERTICES;
$A = new Set<Int>([2,3]);
print rays($A, $O);
$hd = $pc->HASSE_DIAGRAM;
$R = realisation($A, $O);
print old_closure($R, $hd);
print is_face($A,$O,$hd);


application "fan";
$m = matroid::uniform_matroid(3,4);
$t = tropical::matroid_fan<Max>($m);
print $t->COMPACTIFICATION->ALL_FACES;
print $t->COMPACTIFICATION->VERTICES;

application "graph";
$g = complete(4);
application "matroid";
$m = matroid_from_graph($g);
application "tropical";
$t = matroid_fan<Max>($m);
$t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);


application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = divisor( (projective_torus<Max>(2)), rational_fct_from_affine_numerator($f));
application "fan";




