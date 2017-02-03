; # F-cosheaf on tropical hypersurface
application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = divisor((projective_torus<Max>(2)), rational_fct_from_affine_numerator($f));
application "fan";
@betti_usual = ();
@betti_bm = ();
for(my $i=0;$i<2;$i++){
   my $fi = $div->fcosheaf($i);
   my $si=$div->usual_chain_complex($fi);
   my $bmi=$div->borel_moore_complex($fi);
   push @betti_usual, $si->BETTI_NUMBERS;
   push @betti_bm, $bmi->BETTI_NUMBERS;
}
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);

; # F-cosheaf on tropical K3
application "fan";
$k3 = load("../k3.poly");
@numbers = (0..2);
@cosheaves = map{$k3->COSHEAF->[$_]} @numbers;
@usualChainComplexes = map{$k3->usual_chain_complex($_)} @cosheaves;
@bmComplexes = map{$k3->borel_moore_complex($_)} @cosheaves;
@betti_usual = map{$_->BETTI_NUMBERS} @usualChainComplexes; # 30 s
@betti_bm = map{$_->BETTI_NUMBERS} @bmComplexes; # 65 s
print new Matrix(@betti_usual);
print new Matrix(@betti_bm);

$usualChainComplexes[0]->print();
print $usualChainComplexes[0]->BETTI_NUMBERS;



