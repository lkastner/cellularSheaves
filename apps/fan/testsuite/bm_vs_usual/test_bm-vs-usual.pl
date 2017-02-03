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
   my $test = 1;
   for(my $i=0; $i<3; $i++){
      for(my $j=0; $j<3; $j++){
         if($A->elem($i, $j) != $B->elem(2-$i, 2-$j)){
            $test = 0;
         }
      }
   }
   compare_values($tls->name, $test, 1);
}

my $files = `ls .`;
my @tlss = split("\n",$files);
@tlss = grep($_ =~ m/poly/, @tlss);
for my $tls (@tlss){
   compare_bm_usual($tls);
}
