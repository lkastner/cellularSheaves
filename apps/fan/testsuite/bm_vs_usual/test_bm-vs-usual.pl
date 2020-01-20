sub compare_bm_usual{
   my($tlsname) = @_;
   my $tls = load($tlsname);
   my @result1 = ();
   my @result2 = ();
   for(my $i=0;$i<3;$i++){
      my $fi = $tls->fcosheaf($i);
      my $si=$tls->usual_chain_complex($fi);
      my $bmi=$tls->borel_moore_complex($fi);
      push @result1, new Vector<Int>(topaz::betti_numbers($si));
      push @result2, new Vector<Int>(topaz::betti_numbers($bmi));
   }
   my $A = new Matrix<Int>(@result1);
   my $B = new Matrix<Int>(@result2);
   # With the saved matrices we verify that $A->elem($i, $j) == $B->elem(2-$i, 2-$j)
   compare_data($tls->name."_usual", $A);
   compare_data($tls->name."_borelMoore", $B);
}

my $files = `ls .`;
my $c = 0;
my @tlss = split("\n",$files);
@tlss = grep($_ =~ m/poly/, @tlss);
for my $tls (@tlss){
   compare_bm_usual($tls);
   $c++;
   if($c == 10){ last; }
}
