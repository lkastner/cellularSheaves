my $pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
my @result = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wcomplex($i);
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
my $computed = new Matrix(@result);
my $desired = new Matrix([[1,0,0,0],[0,3,0,0],[0,0,3,0],[0,0,0,1]]);
compare_values("cube3",$desired, $computed);
