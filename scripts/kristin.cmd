application "fan";
$f = new PolyhedralFan(RAYS=>[[1,0,0],[1,1,0],[1,0,1],[1,1,1]], MAXIMAL_CONES=>[[0,1,2],[1,2,3]]);
$pc = new PolyhedralComplex($f);
$tau = $pc->HASSE_DIAGRAM->FACES->[3];
print find_max_containing($tau, $pc);
print build_matrix_f($tau, $pc, 1);



$v = new Vector("0 0 0 0 0 1");
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(2,4)->BASES,
VALUATION_ON_BASES=>$v,N_ELEMENTS=>4);
$tls = tropical::linear_space($val_matroid);
$tls->VISUAL;