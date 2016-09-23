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


# Tropical linear spaces

$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
@result = ();
for(my $i=1;$i<=3;$i++){
   my $wi = $tls->wsheaf($i);
   my $si=$tls->usual_chain_complex($wi);
   push @result, $si->BETTI_NUMBERS;
}  
print new Matrix(@result);

