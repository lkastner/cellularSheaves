application "fan";
$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$v = [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
@result1 = ();
@result2 = ();
@result3 = ();
@result4 = ();
for(my $i=0;$i<3;$i++){
   print $i;
#   my $fi = $tls->fcosheaf($i);
   my $wi = $tls->wsheaf($i);
#   my $si=$tls->usual_chain_complex($fi);
#   my $bmi=$tls->borel_moore_complex($fi);
   my $wsi=$tls->usual_chain_complex($wi);
   $wsi->print();
#   my $wbmi=$tls->borel_moore_complex($wi);
#   push @result1, $si->BETTI_NUMBERS;
#   push @result2, $bmi->BETTI_NUMBERS;
   push @result3, $wsi->BETTI_NUMBERS;
#   push @result4, $wbmi->BETTI_NUMBERS;
}
print new Matrix(@result1);
print new Matrix(@result2);
print new Matrix(@result3);
print new Matrix(@result4);



#other vectors: 

#  (0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0) 
#  (0,0,3,1,2,1,0,1,0,2,1,0,2,0,3,1,3,1,0,0) 
#  (0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0) 
#  (1,0,2,0,0,1,0,0,0,0,1,1,1,0,2,0,0,1,0,0) 
#  (1,0,1,0,0,2,0,0,0,0,1,1,1,0,2,0,0,1,0,0) 
#  (3,0,2,0,0,2,0,2,1,2,2,3,2,0,4,0,0,3,0,1)   
#  (4,0,4,0,0,4,0,3,3,3,4,4,4,0,4,0,0,4,0,3)
