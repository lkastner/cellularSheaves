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
