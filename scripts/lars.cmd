

#############################################
application "fan";
$d = 4;
$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube($d))));
$w2 = $pc->wsheaf(2);
for(my $i = 1; $i<$d; $i++){
   my $wi = $pc->wcosheaf($i);
   my $cs = $pc->borel_moore_complex($wi);
   print $i,": ",topaz::betti_numbers($cs),"\n";
}



###########################################################

application "graph";
$g = complete(4);
application "matroid";
$m = matroid_from_graph($g);
print $m->BASES;
application "tropical";
$t = matroid_fan<Max>($m);
application "fan";
print $t->RAYS;
$pcFan = new PolyhedralComplex($t);
$f1 = $pcFan->fcosheaf(1);


###############################################################################

application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = tropical::divisor( (projective_torus<Max>(2)) , rational_fct_from_affine_numerator($f));
application "fan";
print $div->INTERNAL_BASES;
print $div->ORIENTATIONS;
$f1 = $div->fcosheaf(1);


#######################

application "graph";
$g = complete(4);
application "matroid";
$m = matroid_from_graph($g);
application "tropical";
$t = matroid_fan<Max>($m);
$t->VERTICES;
application "fan";
$berg = new PolyhedralComplex($t);


application "tropical";
$f = toTropicalPolynomial("max(0,x+5,y+3, x+y+9)");
$div = divisor( (projective_torus<Max>(2)), rational_fct_from_affine_numerator($f));
application "fan";




