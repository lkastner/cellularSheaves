application "fan";
$A = unit_matrix(4);
$B = new Matrix([[1,1,1,1],[1,0,0,0]]);
print build_matrix($A,$B);
print build_matrix_cpp($A,$B);


## Testing if dualwcomplex and the chain complex for cosheaves work.

application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(4))));
$w1 = $pc->wcomplex(1);
$wd1 = $pc->dualwcomplex(1);
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
