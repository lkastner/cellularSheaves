## For theses examples h^{p,q} = 0 if p!=q. For h^{p,p} we recover 
## the h-vector of the polytope we took the cone over.

application "fan";
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(3))));
@result = ();
for(my $i=0; $i<4; $i++){
   my $w = $pc->wsheaf($i);
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
print new Matrix(@result);

$pc = new PolyhedralComplex(check_fan_objects(new Cone(simplex(5))));
@result = ();
for(my $i=0; $i<6; $i++){
   my $w = $pc->wsheaf($i);
   push @result, $w->CHAIN_COMPLEX->BETTI_NUMBERS;
}
print new Matrix(@result);


