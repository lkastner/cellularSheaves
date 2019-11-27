my $p10 = load("P10_10_11112030203102103020410213253542435.poly");
for(my $i=0; $i<=$p10->DIM; $i++){
   my $w = $p10->wsheaf($i);
   my $s = $p10->usual_cochain_complex($w);
   my $array = convert_cc($s);
   compare_data("P10_UCCC_$i", $array);
} 

for(my $i=0;$i<3;$i++){
   my $fi = $p10->fcosheaf($i);
   my $si=$p10->usual_chain_complex($fi);
   my $bmi=$p10->borel_moore_complex($fi);
   my $array_s = convert_cc($si);
   my $array_bm = convert_cc($bmi);
   compare_data("P10_UCC_$i", $array_s);
   compare_data("P10_BMC_$i", $array_bm);
}


my $k3 = load("k3.pcom");
for(my $i=0; $i<3; $i++){
   my $f = $k3->compact_fcosheaf($i);
   my $d = new topaz::ChainComplex<Matrix<Rational>>(build_full_chain($k3->COMPACTIFICATION, $k3->COMPACTIFICATION->ORIENTATIONS, $f, false));
   my $array = convert_cc($d);
   compare_data("k3_$i", $array);
}

