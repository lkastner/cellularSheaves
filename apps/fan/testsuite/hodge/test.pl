my $ap = application "topaz";

my $k3 = load("k3.pcom");
my @rows_k3;
for(my $i=0; $i<3; $i++){
   my $f = $k3->compact_fcosheaf($i);
   my $d = build_full_chain($k3->COMPACTIFICATION, $k3->COMPACTIFICATION->ORIENTATIONS, $f, false);
   push @rows_k3, new Vector<Int>(topaz::betti_numbers(new topaz::ChainComplex<Matrix<Rational>>($d)));
}
my $computed = new Matrix(@rows_k3);
my $desired = new Matrix([[1,0,1],[0,20,0],[1,0,1]]);
compare_values("k3", $computed, $desired);


my $cube2 = new PolyhedralComplex(check_fan_objects(cube(2)));
my @rows_cube2;
for(my $i=0; $i<3; $i++){
   my $f = $cube2->compact_fcosheaf($i);
   my $d = build_full_chain($cube2->COMPACTIFICATION, $cube2->COMPACTIFICATION->ORIENTATIONS, $f, false);
   push @rows_cube2, new Vector<Int>(topaz::betti_numbers(new topaz::ChainComplex<Matrix<Rational>>($d)));
}
$computed = new Matrix(@rows_cube2);
$desired = new Matrix([[1,0,0],[2,0,0],[1,0,0]]);
compare_values("cube2", $computed, $desired);

my $cube3 = new PolyhedralComplex(check_fan_objects(cube(3)));
my @rows_cube3;
for(my $i=0; $i<4; $i++){
   my $f = $cube3->compact_fcosheaf($i);
   my $d = build_full_chain($cube3->COMPACTIFICATION, $cube3->COMPACTIFICATION->ORIENTATIONS, $f, false);
   push @rows_cube3, new Vector<Int>(topaz::betti_numbers(new topaz::ChainComplex<Matrix<Rational>>($d)));
}
$computed = new Matrix(@rows_cube3);
$desired = new Matrix([[1,0,0,0],[3,0,0,0],[3,0,0,0],[1,0,0,0]]);
compare_values("cube3", $computed, $desired);

my $cubicSurface = load("cubicSurface.pcom");
my @rows_cubicSurface;
for(my $i=0; $i<3; $i++){
   my $f = $cubicSurface->compact_fcosheaf($i);
   my $d = build_full_chain($cubicSurface->COMPACTIFICATION, $cubicSurface->COMPACTIFICATION->ORIENTATIONS, $f, false);
   push @rows_cubicSurface, new Vector<Int>(topaz::betti_numbers(new topaz::ChainComplex<Matrix<Rational>>($d)));
}
$computed = new Matrix(@rows_cubicSurface);
$desired = new Matrix([[1,0,0],[0,7,0],[0,0,1]]);
compare_values("cubicSurface", $computed, $desired);
