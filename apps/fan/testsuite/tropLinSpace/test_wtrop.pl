my $tls = load("00100010000000101000.poly");
my $computed = new Matrix([[0,0,0]]);
for(my $i=0;$i<3;$i++){
   print $i;
   my $wi = $tls->wsheaf($i);
   my $wsi=$tls->usual_chain_complex($wi);
   $computed->elem(0,$i) = $wsi->IS_WELLDEFINED;
}
my $desired = new Matrix([[1,1,1]]);
compare_values("WTrop",$desired,$computed);
