// center: [0, tan(alpha)]
// radius: 2 - sec(alpha)
// alpha: [0, pi/4]
include <params.scad>
// include <icosphere.scad>
include <egg.scad>

module support_ball_holes() {
  for(i = [0:2])
    rotate([0, 0, 120 * i])
      rotate([-45, 0, 0])
        translate([0, 0, -hole_r])
          scale(1.25)
            // import("icosphere-3.stl");
            cylinder(h = 2.2, r = 1, center = true, $fn = 64);
}

module shell_base() {
  difference() {
    rotate([egg_tilt, 0, 0])
      translate(-center)
        translate([0, 0, btm_h])
          rotate([-90, 0, 0])
            scale(egg_scale)
              // import("egg-42-4.stl");
              egg(egg_btm_alpha, egg_top_alpha, 4);
    scale(hole_r)
      import("icosphere-4.stl");
    translate([0, 0, -200 - center[2]])
      cube(400, center = true);
  }
}

module shell_3d_minkowski() {
  difference() {
    minkowski() {
      difference() {
        rotate([egg_tilt, 0, 0])
          translate(-center)
            translate([0, 0, btm_h])
              rotate([-90, 0, 0])
                scale(egg_mkw_scale)
                  import("egg-42-2.stl");
        scale(hole_mkw_r)
          import("icosphere-2.stl");
      }
      scale(mkw_r)
        import("icosphere-1.stl");
    }
    translate([0, 0, -200 - center[2]])
      cube(400, center = true);
  }
}

tilt = -24;

module egg_base() {
  rotate([tilt, 0, 0])
    translate(-center)
      rotate([-90, 0, 0])
        scale(egg_mkw_scale)
          import("egg-42-4.stl");
}

dazim = 12;
delev = 12;

r = 300;
th = delev / 2;
x = r * sin(th);
y = r * cos(th);
module elev_triangle() {
  pnts = [
    [0, 0], 
    [+x, y], 
    [-x, y]
  ];
  polygon(pnts);
}

module azim_slice(azim) {
  intersection() {
    offset(mkw_r)
      difference() {
        projection(cut = true)
          rotate([0, 90, 0])
            rotate([0, 0, azim])
              children();
        circle(r = hole_mkw_r);
      }
  }
}

module elev_section(elev, thick = 1) {
  rotate([0, 0, -elev])
    translate([0, 0, -thick / 2])
      linear_extrude(thick, scale = 1)
        intersection() {
          rotate([0, 0, elev])
            children();
          elev_triangle();
        }
}

module half_track_egg() {
  rotate([-tilt, 0, 0])
    union()
      for(azim0 = [dazim / 2:dazim:180])
        union()
          for(elev = [-90 + delev / 2:delev / 2:90])
            hull()
              for(sgn = [-1, +1]) {
                azim = azim0 + sgn * dazim / 2;
                rotate([0, 0, -azim])
                  rotate([0, -90, 0])
                    elev_section(elev)
                      azim_slice(azim)
                        egg_base();
              }
}

module shell_2d_offset() {
  difference() {
    rotate([egg_tilt, 0, 0])
      union() {
        half_track_egg();
        mirror([1, 0, 0])
          half_track_egg();
      }
    translate([0, 0, -200 - center[2]])
      cube(400, center = true);
  }
}

module shell_vtk(path) {
  difference() {
    translate([0, 0, btm_h])
      // rotate([0, 0, 180])
      rotate([0, 180, 0])
        rotate([-90, 0, 0])
          import(path);
    translate([0, 0, -200])
      cube(400, center = true);
  }
}

// shell
if (true) {
// shell_base();
// shell_3d_minkowski();
// shell_2d_offset();
// shell_vtk("../../surface-mkw=3_mc.stl");
// shell_vtk("../../surface-mkw=3_fm.stl");
}
difference() {
  shell_base();
  // shell_3d_minkowski();
  // shell_2d_offset();
  // shell_vtk("../../surface-mkw=4_fm.stl");
  support_ball_holes();
}

// ball
if (false)
  color("white")
    translate(center)
      scale(ball_r)
        import("icosphere-4.stl");