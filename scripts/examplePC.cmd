application "fan";
$f = new PolyhedralFan(RAYS=>[[1,0,0],[1,1,0],[1,0,1],[1,1,1]], MAXIMAL_CONES=>[[0,1,2],[1,2,3]]);
$pc = new PolyhedralComplex($f);
print keys %{$pc->CHOSEN_BASES};

$s = new Set<Int>(1,2,3);
print $pc->CHOSEN_BASES->{$s};

$pc->BLOCKS;
$h = new Set<Int>([0]);
print $pc->CHOSEN_BASES->{$h}->type->full_name;


$t = new Set<Int>(2,3);
$pst = new Pair<Set<Int>, Set<Int> >($s, $t);
print $pc->BLOCKS->{$pst};
print $pc->CHOSEN_BASES->{$s};
print $pc->CHOSEN_BASES->{$t};

$pc->CHAIN_COMPLEX;

#############################################

$pc = new PolyhedralComplex(check_fan_objects(new Cone(cube(4))));
$pc->CHAIN_COMPLEX;


