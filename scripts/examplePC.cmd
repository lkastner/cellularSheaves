application "fan";
$f = new PolyhedralFan(RAYS=>[[1,0,0],[1,1,0],[1,0,1],[1,1,1]], MAXIMAL_CONES=>[[0,1,2],[1,2,3]]);
$pc = new PolyhedralComplex($f);
print keys %{$pc->CHOSEN_BASES};

$s = new Set<Int>(1,2,3);
print $pc->CHOSEN_BASES->{$s};

$pc->BLOCKS;
$h = new Set<Int>([0]);
print $pc->CHOSEN_BASES->{$h}->type->full_name;


$t = new Set<Int>(2,3);
$pst = new Pair<Set<Int>, Set<Int> >($s, $t);
print $pc->BLOCKS->{$pst};
print $pc->CHOSEN_BASES->{$s};
print $pc->CHOSEN_BASES->{$t};

$pc->CHAIN_COMPLEX;

#############################################
application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(4))));
$pc->CHAIN_COMPLEX;
print $pc->CHAIN_COMPLEX->[2]*$pc->CHAIN_COMPLEX->[3];
print $pc->CHAIN_COMPLEX->[3]*$pc->CHAIN_COMPLEX->[4];

$A = $pc->CHAIN_COMPLEX->[2];
$B = $pc->CHAIN_COMPLEX->[3];
print compute_homology($B,$A);
print compute_homology_dimension($B,$A);


$A = new Matrix([[2, 1/2, 0],[1,0,0]]);
print dense(to_lattice_basis($A));


application "fan";
$A = new Matrix([[1,2,3,4],[3,4,5,6],[7,8,9,0]]);
print wedge_matrix($A, 2);

application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
$pc->CHAIN_COMPLEX;
$c2 = $pc->CHAIN_COMPLEX->wedge_complex(2);
$w2 = $pc->wcomplex(2);

application "fan";
$d = 4;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$pc->CHAIN_COMPLEX;
$w2 = $pc->wcomplex(2);
$s = new Set<Int>(12,13,14,15);
$t = new Set<Int>(14,15);
$p = new Pair<Set<Int>, Set<Int> >($s, $t);
print $pc->BLOCKS->{$p};
print $w2->BLOCKS->{$p};
print $w2->CHAIN_COMPLEX->IS_WELLDEFINED;
print wedge_matrix($pc->BLOCKS->{$p},2)->cols;


for(my $i = 1; $i<$d; $i++){
   my $wi = $pc->wcomplex($i);
   print $i,": ",$wi->CHAIN_COMPLEX->BETTI_NUMBERS,"\n";
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
$blocks = $pc->BLOCKS;
print assemble_matrix($sigmas, $taus, $blocks, $orientations);


###########################################################

#Testing f blocks


print build_generator_matrix_f($tau, $pc, 1); 