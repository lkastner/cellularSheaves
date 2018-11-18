my $berg = load("completegraph4.poly");
my @result1 = ();
my @result2 = ();
for(my $i=0; $i<3; $i++){
   my $f = $berg->fcosheaf($i);
   my $s = $berg->usual_chain_complex($f);
   my $bm = $berg->borel_moore_complex($f);
   push @result1, $s->BETTI_NUMBERS;
   push @result2, $bm->BETTI_NUMBERS;
}
my $computed1 = new Matrix<Int>(@result1);
my $computed2 = new Matrix<Int>(@result2);
my $desired1 = new Matrix<Int>([1,0,0],[5,0,0],[6,0,0]);
my $desired2 = new Matrix<Int>([0,0,6],[0,0,5],[0,0,1]);
compare_values("completegraph4usual",$computed1,$desired1);
compare_values("completegraph4bm",$computed2,$desired2);
