my $pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
my @result = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wcosheaf($i);
   my $cs = $pc->borel_moore_complex($w);
   push @result, $cs->BETTI_NUMBERS;
}
my $computed = new Matrix<Int>(@result);
my $desired = new Matrix<Int>([[1,0,0,0],[0,3,0,0],[0,0,3,0],[0,0,0,1]]);
compare_values("cube3",$desired, $computed);
