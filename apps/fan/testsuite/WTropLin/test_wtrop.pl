application "fan";
$v = [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
$computed = new Matrix([[0,0,0]]);
for(my $i=0;$i<3;$i++){
   print $i;
   my $wi = $tls->wsheaf($i);
   my $wsi=$tls->usual_chain_complex($wi);
   $computed->elem(0,$i) = $wsi->IS_WELLDEFINED;
}
$desired = new Matrix([[1,1,1]]);
compare_values("WTrop",$desired,$computed);
