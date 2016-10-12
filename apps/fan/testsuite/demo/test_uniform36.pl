my $berg = load("uniform36.poly");
my @result1 = ();
my @result2 = ();
for(my $i=0; $i<6; $i++){
   my $f = $berg->fcosheaf($i);
   my $s = $berg->usual_chain_complex($f);
   my $bm = $berg->borel_moore_complex($f);
   push @result1, $s->BETTI_NUMBERS;
   push @result2, $bm->BETTI_NUMBERS;
}
my $computed1 = new Matrix(@result1);
my $computed2 = new Matrix(@result2);
my $desired1 = new Matrix([1,0,0],[5,0,0],[10,0,0],[0,0,0],[0,0,0],[0,0,0]);
my $desired2 = new Matrix([0,0,10],[0,0,5],[0,0,1],[0,0,0],[0,0,0],[0,0,0]);
compare_values("uniform36usual",$computed1,$desired1);
compare_values("uniform36bm",$computed2,$desired2);
