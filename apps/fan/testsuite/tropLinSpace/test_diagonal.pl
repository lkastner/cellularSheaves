sub is_diagonal_matrix{
   my ($A)= @_;
   for (my $r = 0; $r<$A->rows;$r++){
      for (my $c = 0; $c<$A->cols;$c++){
         if($r!=$c && $A->elem($r,$c)!=0){
            return 0;
         }
      }
   }
   return 1;
}

sub test_tls{
   my ($tls) = @_;
   my @result = ();
   for(my $i=0;$i<=$tls->DIM;$i++){
      my $wi = $tls->wsheaf($i);
      my $si=$tls->usual_chain_complex($wi);
      push @result, $si->BETTI_NUMBERS;
   }
   my $A = new Matrix(@result);
   my $test = is_diagonal_matrix($A);
   compare_values($tls->name,$test,1);
}

## build the tropical space in polymake and save it.
my @tlss = ("00100010000000101000.poly",  "10200100001110200100.poly",  
"00312101021020313100.poly",  "30200202122320400301.poly",  
"00312101022030412200.poly",  "40400403334440400403.poly",
"10100200001110200100.poly");
for my $tls (@tlss){
   my $tlsReal = load($tls);
   test_tls($tlsReal);
}
