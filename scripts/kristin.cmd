application "fan";
$f = new PolyhedralFan(RAYS=>[[1,0,0],[1,1,0],[1,0,1],[1,1,1]], MAXIMAL_CONES=>[[0,1,2],[1,2,3]]);
$pc = new PolyhedralComplex($f);
$tau = $pc->HASSE_DIAGRAM->FACES->[3];
$sigma = $pc->HASSE_DIAGRAM->FACES->[2];
print find_max_containing($tau, $pc);
print build_matrix_f($tau, $pc, 1);
print build_blocks_f($pc, 1);
$f2 = $pc->fcomplex(2);
print $f2->CHAIN_COMPLEX->DIFFERENTIALS;

$key = new Pair<Set<Int>, Set<Int> >($sigma, $tau);
print $key;
print build_blocks_f($pc, 1)->{$key};


$v = new Vector("0 0 0 0 0 1");
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(2,4)->BASES,
VALUATION_ON_BASES=>$v,N_ELEMENTS=>4);
$tls = tropical::linear_space($val_matroid);
$tls->VISUAL;


application "fan";
$poly = new Polytope(POINTS=>[[0, 1, 0], [1,0,0], [0, 0, 1]]);
$pc = new PolyhedralComplex(check_fan_objects(new Cone($poly)));
print $pc->HASSE_DIAGRAM->FACES;
print $pc->BOUNDED_FACES;
print $pc->UNBOUNDED_FACES;
print $pc->FAR_FACES;
print $pc->NON_FAR_FACES;
$w4 = $pc->wcomplex(4);
$w3 = $pc->wcomplex(3);
$w2 = $pc->wcomplex(2);
$w1 = $pc->wcomplex(1);
$boundedChain = build_chain_complex($w4->BLOCKS, $pc->BOUNDED_FACES, $pc->ORIENTATIONS);
$boundedChain = build_chain_complex($w3->BLOCKS, $pc->BOUNDED_FACES, $pc->ORIENTATIONS);
$boundedChain = build_chain_complex($w2->BLOCKS, $pc->BOUNDED_FACES, $pc->ORIENTATIONS);
$boundedChain = build_chain_complex($w1->BLOCKS, $pc->BOUNDED_FACES, $pc->ORIENTATIONS);


###printing all differential dimensions.

map{print $_->rows, " ", $_->cols, "\n"}@{$bm2->DIFFERENTIALS};


#####################
application "matroid";
$m = uniform_matroid(2,4); 
application "tropical";
$mFan = matroid_fan<Max>($m);
application "fan";
$pcFan = new PolyhedralComplex($mFan);
$f1 = $pcFan->fcomplex(1);

application "tropical";
$v = new matroid::ValuatedMatroid<Min>(N_ELEMENTS=>5, BASES=>[[0,1,2],[0,1,3],[0,1,4],[0,2,3],[0,2,4],[0,3,4], [1,2,3], [1,2,4],[1,3,4],[2,3,4]], VALUATION_ON_BASES=>[1,0,3,2,6,1,0,3,1,0]);
$l = tropical::linear_space($v);


$v3 = new matroid::ValuatedMatroid<Min>(N_ELEMENTS=>6, BASES=>[[0,1,2],[0,1,3],[0,1,4],[0,1,5],[0,2,3],[0,2,4], [0,2,5], [0,3,4],[0,3,5], [0,4,5], [1,2,3], [1,2,4], [1,2,5],[1,3,4],[1,3,5], [1,4,5], [2,3,4], [2,3,5], [2,4,5],[3,4,5]], VALUATION_ON_BASES=>[0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0]);
$l3 = tropical::linear_space($v3);
print $l3->F_VECTOR;
$p = 1;
$boundedChain = build_chain_complex($l->$wcomplex($p)->BLOCKS, $l->BOUNDED_FACES, $l->ORIENTATIONS);

$v2 = new matroid::ValuatedMatroid<Min>(N_ELEMENTS=>5, BASES=>[[3,4],[2,4],[2,3],[1,4],[1,3],[1,2], [0,4], [0,3],[0,2],[0,1]], VALUATION_ON_BASES=>[1,0,3,2,6,1,0,3,1,0]);


application "tropical";
$f = toTropicalPolynomial("max(0,x,y,z, 2*x +1, 2*y+5, x*y-8)");
$div = divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
$pc = new PolyhedralComplex($div);
print $pc->VERTICES;
$f1= $pc->fcomplex(1);




application "tropical";
$f = toTropicalPolynomial("max(0,x,y,z)");
$div = divisor( (projective_torus<Max>(3)) , rational_fct_from_affine_numerator($f));
application "fan";
$f1 = $div->fcomplex(1);
$f2  = $div->fcomplex(2);
$f0 = $div->fcomplex(0);
$f3 = $div->fcomplex(3);
$bm1  = $div->borel_moore_complex($f1); 
$bm2  = $div->borel_moore_complex($f2); 
$bm0  = $div->borel_moore_complex($f0); 
$bm3  = $div->borel_moore_complex($f3); 
print $bm1->IS_WELLDEFINED;
print $bm2->IS_WELLDEFINED;
print $bm3->IS_WELLDEFINED;




$pc = new PolyhedralComplex($div);
print $pc->VERTICES;
$f1= $pc->fcomplex(1);



application "tropical";
$f = toTropicalPolynomial("max(0,x,y)");
$div = divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
$pc = new PolyhedralComplex($div);
print $pc->VERTICES;
$f1= $pc->fcomplex(1);


####testing constant_sheaf

application "fan";
$d = 2;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$const = $pc->constant_sheaf;
$bmcomp = $pc->borel_moore_complex($const);
print $bmcomp->DIFFERENTIALS;
print $bmcomp->BETTI_NUMBERS;