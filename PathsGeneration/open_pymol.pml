# Set color space to CMYK
space cmyk

# Apply performance settings
set line_smooth, on
set depth_cue, on
set specular, 1.0
set surface_quality, 1
set cartoon_sampling, 14
set ribbon_sampling, 10
set transparency_mode, 2
set use_shaders, on
set cartoon_use_shader, on
set cgo_use_shader, on
set dash_use_shader, on
set dot_use_shader, on
set line_use_shader, on
set mesh_use_shader, on
set nb_spheres_use_shader, 1
set nonbonded_use_shader, on
set ribbon_use_shader, on
set sphere_use_shader, on
set stick_use_shader, on
set surface_use_shader, on
set render_as_cylinders, on
set alignment_as_cylinders, on
set cartoon_nucleic_acid_as_cylinders, 1
set dash_as_cylinders, on
set line_as_cylinders, on
set mesh_as_cylinders, on
set nonbonded_as_cylinders, on
set ribbon_as_cylinders, on
set stick_as_cylinders, on
set dot_as_spheres, on
set stick_ball, off
set sphere_mode, 9
set nb_spheres_quality, 3

# Apply rebuild settings
rebuild
set specular, 0.0
set stick_ball, on
set stick_ball_ratio, 1.5
set stick_radius, 0.2
show_as sticks

# Color
util.cbag

# Mouse
set mouse_selection_mode, 0

# Labels
set label_size, 24
set label_font_id, 7
# Label all atoms except hydrogens
label not elem H, name