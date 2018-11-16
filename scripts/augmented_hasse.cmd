

#######################
# C++ Compactification
application "fan";
$pc = new PolyhedralComplex(POINTS=>[[1,0,0],[0,1,0],[0,0,1]],INPUT_POLYTOPES=>[[0,1,2]]);
tropcomp($pc);
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[1,1,0],[1,0,1],[0,1,1]], INPUT_CONES=>[[0,1,3],[1,2,3]]);
tropcomp($pc);

application "fan";
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0,0],[0,1,0,0],[0,1,1,0],[0,1,1,1],[0,1,0,1]], INPUT_CONES=>[[0,1,2,3,4]]);
$hd = tropcomp($pc);

application "fan";
# Positive 2-dim orthant
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[0,1,0],[0,0,1]], INPUT_CONES=>[[0,1,2]]);
$hd = tropcomp($pc);
# Positive 2-dim orthant without interior
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[0,1,0],[0,0,1]], INPUT_CONES=>[[0,1],[0,2]]);
$hd = tropcomp($pc);
# Thickened x-axis
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[1,0,1],[0,1,0]], INPUT_CONES=>[[0,1,2]]);
$hd = tropcomp($pc);
# Thickened x-axis without interior
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[1,0,1],[0,1,0]], INPUT_CONES=>[[0,1],[0,2],[1,2]]);
$hd = tropcomp($pc);

application "fan";
# Positive 2-dim orthant
$pc = new PolyhedralComplex(INPUT_RAYS=>[[1,0,0],[0,1,0],[0,0,1]], INPUT_CONES=>[[0,1,2]]);
$hd = tropcomp($pc);
$f1 = $pc->fcosheaf(1);
$cp = $pc->usual_chain_complex($f1);

application "fan";
$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);
$tls = tropical::linear_space($val_matroid);
$f1 = $tls->fcosheaf(1);
$cp = $tls->usual_chain_complex($f1);

