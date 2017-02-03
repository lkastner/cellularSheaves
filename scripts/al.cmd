application "fan";
$A = unit_matrix(4);
$B = new Matrix([[1,1,1,1],[1,0,0,0]]);
print build_matrix($A,$B);
print build_matrix_cpp($A,$B);


## Testing if dualwcomplex and the chain complex for cosheaves work.

application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(4))));
$w1 = $pc->wsheaf(1);
$wd1 = $pc->wcosheaf(1);
print $wd1->CHAIN_COMPLEX->IS_WELLDEFINED;
for(my $i = 0; $i<=7; $i++){
   print $i;
   print " : ";
   print $wd1->CHAIN_COMPLEX->DIFFERENTIALS->[$i]->rows;
   print ", ";
   print $wd1->CHAIN_COMPLEX->DIFFERENTIALS->[$i]->cols;
   print ". \n";
}
print $wd1->CHAIN_COMPLEX->BETTI_NUMBERS;
for(my $i = 0; $i<=7; $i++){
   print $i;
   print " : ";
   print $w1->CHAIN_COMPLEX->DIFFERENTIALS->[$i]->rows;
   print ", ";
   print $w1->CHAIN_COMPLEX->DIFFERENTIALS->[$i]->cols;
   print ". \n";
}
print $w1->CHAIN_COMPLEX->BETTI_NUMBERS;
$w2 = $pc->wsheaf(2);
$wd2 = $pc->wcosheaf(2);
print $w2->CHAIN_COMPLEX->BETTI_NUMBERS;
print $wd2->CHAIN_COMPLEX->BETTI_NUMBERS;


###############################################################
############### Problems with zero maps again #################
###############################################################
application "fan";
$v = [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
$w0 = $tls->wsheaf(0);
$w1 = $tls->wsheaf(1);
$w2 = $tls->wsheaf(2);
$ws0 = $tls->usual_chain_complex($w0);
$ws1 = $tls->usual_chain_complex($w1);
$ws2 = $tls->usual_chain_complex($w2);
$ws0->print_debug();
$ws1->print_debug();
$ws2->print_debug();


application "fan";
$v = [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
$w1 = $tls->wsheaf(1);
$ws1 = $tls->usual_chain_complex($w1);
$ws1->print_debug();


###############################################################
############### Checking print method for cochain complex #####
###############################################################
application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(4))));
$w1 = $pc->wsheaf(1);
$wd1 = $pc->wcosheaf(1);
$s1 = $pc->usual_cochain_complex($w1);
$sd1 = $pc->usual_chain_complex($wd1);
$sd1->print();
$s1->print();
