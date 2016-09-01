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
   my $w = $pc->wcomplex($i);
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
print new Matrix(@result);

application "fan";
$d = 4;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$w2 = $pc->wcomplex(2);
$s = new Set<Int>(12,13,14,15);
$t = new Set<Int>(14,15);
$p = new Pair<Set<Int>, Set<Int> >($s, $t);
print $w2->BLOCKS->{$p};
print wedge_matrix($pc->SIMPLE_BLOCKS->{$p},2)->cols;


for(my $i = 1; $i<$d; $i++){
   my $wi = $pc->wcomplex($i);
   print $i,": ",$wi->BETTI_NUMBERS,"\n";
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
$f1 = $pcFan->fcomplex(1);
$boundedChain = build_chain_complex($f1->BLOCKS, $pcFan->BOUNDED_FACES, $pcFan->ORIENTATIONS);



###############################################################################


application "tropical"; 
$f = toTropicalPolynomial("max(0,x,y,z)"); 
$div = divisor( (projective_torus<Max>(3)) , rational_fct_from_affine_numerator($f)); 
application "fan"; 
$f1 = $div->fcomplex(1); 
$bm1 = $div->borel_moore_complex($f1); 
$bm1->print();

$f2 = $div->fcomplex(2);
$f0 = $div->fcomplex(0); 
$f3 = $div->fcomplex(3); 
$bm1 = $div->borel_moore_complex($f1); 
$bm2 = $div->borel_moore_complex($f2); 
$bm0 = $div->borel_moore_complex($f0); 
$bm3 = $div->borel_moore_complex($f3); 
$bm1->print();
print $bm1->IS_WELLDEFINED;

###############################################################################

application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = tropical::divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
print $div->ORIENTATIONS;
$f1 = $div->fcomplex(1);



