
###########################################################
application "fan";
$pc = new PolyhedralComplex(VERTICES=>[[1,0,0],[0,-1,0],[0,0,-1],[0,1,1]], MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3]]);
print $pc->SIMPLE_BLOCKS;
$w = $pc->wcomplex(1);
$w->CHAIN_COMPLEX->print();
print $w->CHAIN_COMPLEX->DIFFERENTIALS;
print $w->CHAIN_COMPLEX->BETTI_NUMBERS;
$f1 = $pc->fcomplex(1);
$boundedChain = build_chain_complex($f1->BLOCKS, $pc->BOUNDED_FACES, $pc->ORIENTATIONS);

