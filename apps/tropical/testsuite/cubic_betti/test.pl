foreach my $file (glob("*.pcom")){
   my ($name) = $file =~ m/(.*)\.pcom/;
   my $trop = load($file);
   my $comp = $trop->COMPACTIFICATION;
   my $pw = $trop->PATCHWORK;
   my $cosheaf = $pw->sign_cosheaf();
   my $chain = fan::build_full_chain($comp, $comp->ORIENTATIONS, $cosheaf, true);
   my $betti = topaz::betti_numbers<GF2>($chain);
   compare_data($name, $betti);
}
