load example.pdb, prot
load example_coms.pdb, coms

hide everything, prot
show cartoon, prot

hide everything, coms
show spheres, coms
color black, coms

# Draw distances
distance dist1, coms///1/CA, coms///3/CA
distance dist2, coms///2/CA, coms///3/CA
distance dist3, coms///1/CA, coms///2/CA

# Show angle
angle ang1, coms///1/CA, coms///3/CA, coms///2/CA

set dash_color, black
set dash_width, 2
set sphere_scale, 0.6, coms

bg_color white
