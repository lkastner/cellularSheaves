
use LWP::Simple;

# open(FILE, "saved.html") or die "$!";
# local $/;
# $document = <FILE>;
# close FILE;

$document = "LWP::Simple".get("http://www.uni-math.gwdg.de/jensen/Research/G3_7/grassmann3_7.html");


@tables = $document =~ m/\<table.*?\/table\>/gs;
@rows = $tables[0] =~ m/\<tr.*?\/tr\>/gs;
@rows = grep($_ =~ m/\([-?\d+,]{50,}\)/gs, @rows);
$i = 0;
$currentName = "";
for my $row (@rows) {
   if($row =~ m/(P\d+)/g){
      my $name = $1;
      $currentName = $name;
   }
   my($pluecker) = $row =~ m/\([-?\d+,]{50,}\)/gs;
   print $i,": Name: ",$currentName," Pluecker: ",$pluecker,"\n";
   my @v = $pluecker =~ m/(-?\d+)/g;
   my $v = new Vector(\@v);
   my $saveName = $currentName."_".$i."_".join("",@v);
   print $saveName,"\n";
   $i++;
   my $val_matroid = new matroid::ValuatedMatroid<Min>(BASES=>matroid::uniform_matroid(3,7)->BASES,VALUATION_ON_BASES=>$v,N_ELEMENTS=>7);
   print matroid::check_valuated_basis_axioms($val_matroid->BASES,$val_matroid->VALUATION_ON_BASES,verbose=>1);
   my $tls = tropical::linear_space($val_matroid);
   $tls->name = $saveName;
   save($tls, $saveName.".poly");
}
