my $k3 = load("k3.pcom");
compare_data("k3", hodge_numbers($k3));

my $cube2 = new PolyhedralComplex(check_fan_objects(cube(2)));
compare_data("cube2", hodge_numbers($cube2));

my $cube3 = new PolyhedralComplex(check_fan_objects(cube(3)));
compare_data("cube3", hodge_numbers($cube3));

my $cubicSurface = load("cubicSurface.pcom");
compare_data("cubicSurface", hodge_numbers($cubicSurface));
