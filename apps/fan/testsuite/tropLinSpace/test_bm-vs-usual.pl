sub compare_bm_usual{
   my($tlsname) = @_;
   my $tls = load($tlsname);
   my @result1 = ();
   my @result2 = ();
   for(my $i=0;$i<3;$i++){
      my $fi = $tls->fcosheaf($i);
      my $si=$tls->usual_chain_complex($fi);
      my $bmi=$tls->borel_moore_complex($fi);
      push @result1, $si->BETTI_NUMBERS;
      push @result2, $bmi->BETTI_NUMBERS;
   }
   my $A = new Matrix(@result1);
   my $B = new Matrix(@result2);
   for(my $i=0; $i<3; $i++){
      for(my $j=0; $j<3; $j++){
         compare_values($tls->name."(".$i.",".$j.")", $A->elem($i, $j), $B->elem(2-$i, 2-$j));
      }
   }
}

my @tlss = ("00100010000000101000.poly",  "10200100001110200100.poly",  
"00312101021020313100.poly",  "30200202122320400301.poly",  
"00312101022030412200.poly",  "40400403334440400403.poly",
"10100200001110200100.poly");
for my $tls (@tlss){
   compare_bm_usual($tls);
}
