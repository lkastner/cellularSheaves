
###########################################################
application "fan";
$pc = new PolyhedralComplex(VERTICES=>[[1,0,0],[0,-1,0],[0,0,-1],[0,1,1]], MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3]]);
print $pc->SIMPLE_BLOCKS;
@result = ();
for(my $i=0; $i<4; $i++){
   $w = $pc->wsheaf($i);
   print $w->CHAIN_COMPLEX->BETTI_NUMBERS;
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
print new Matrix(@result);
$f1 = $pc->fcosheaf(1);
$boundedChain = build_chain_complex($f1->BLOCKS, $pc->BOUNDED_FACES, $pc->ORIENTATIONS);

