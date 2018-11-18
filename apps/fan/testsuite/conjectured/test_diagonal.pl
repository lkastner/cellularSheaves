sub test_tls{
   my ($tls) = @_;
   my @result = ();
   for(my $i=0;$i<=$tls->DIM;$i++){
      my $wi = $tls->wsheaf($i);
      my $si=$tls->usual_cochain_complex($wi);
      push @result, $si->BETTI_NUMBERS;
   }
   my $A = new Matrix<Int>(@result);
   compare_data($tls->name, $A);
}

## build the tropical space in polymake and save it.
my $files = `ls .`;
my @tlss = split("\n",$files);
@tlss = grep($_ =~ m/poly/, @tlss);
for my $tls (@tlss){
   my $tlsReal = load($tls);
   test_tls($tlsReal);
}
