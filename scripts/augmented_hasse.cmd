

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
