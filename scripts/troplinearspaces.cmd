$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(2,4)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>4);


$v = [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0];


$v = [0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0];
$val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,6)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>6);

print matroid::check_valuated_basis_axioms($val_matroid->BASES,$val_matroid->VALUATION_ON_BASES,verbose=>1);

$tls = tropical::linear_space($val_matroid);


#other vectors: 

#	(0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0)	
#	(0,0,3,1,2,1,0,1,0,2,1,0,2,0,3,1,3,1,0,0)	
#	(0,0,3,1,2,1,0,1,0,2,2,0,3,0,4,1,2,2,0,0)	
#	(1,0,2,0,0,1,0,0,0,0,1,1,1,0,2,0,0,1,0,0)	
#	(1,0,1,0,0,2,0,0,0,0,1,1,1,0,2,0,0,1,0,0)	
#(3,0,2,0,0,2,0,2,1,2,2,3,2,0,4,0,0,3,0,1)	
#	(4,0,4,0,0,4,0,3,3,3,4,4,4,0,4,0,0,4,0,3)	