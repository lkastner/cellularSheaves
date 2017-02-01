; # Vanishing of W-sheaf cohomology
@wbetti_usual = ();
@wbetti_cs = ();
for(my $i=0;$i<3;$i++){
   my $wi = $tls->wsheaf($i);
   my $wsi=$tls->usual_cochain_complex($wi);
   my $wcsi=$tls->compact_support_complex($wi);
   push @wbetti_usual, $wsi->BETTI_NUMBERS;
   push @wbetti_cs, $wcsi->BETTI_NUMBERS;
}
print new Matrix(@wbetti_usual);
print new Matrix(@wbetti_cs);

load_commands("tropical.cmd");
