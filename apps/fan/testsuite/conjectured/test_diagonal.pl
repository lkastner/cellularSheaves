sub test_tls{
   my ($tls) = @_;
   my @result = ();
   for(my $i=0;$i<=$tls->DIM;$i++){
      my $wi = $tls->wsheaf($i);
      my $si=$tls->usual_cochain_complex($wi);
      push @result, new Vector<Int>(topaz::betti_numbers($si));
   }
   my $A = new Matrix<Int>(@result);
   compare_data($tls->name, $A);
}

## build the tropical space in polymake and save it.
my $c = 0;
for my $tls (glob("*.poly")){
   my $tlsReal = load($tls);
   test_tls($tlsReal);
   $c++;
   if($c == 10){ last; }
}
