my $tls = load("00312101022030412200.poly");
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
my $C = new Matrix([[1,0,0],[5,0,0],[10,0,0]]);
my $D = new Matrix([[0,0,10],[0,0,5],[0,0,1]]);
compare_values($tls->name."Usual",$A, $C);
compare_values($tls->name."BM",$B, $D);

