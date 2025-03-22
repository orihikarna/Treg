// center: [0, tan(alpha)]
// radius: 2 - sec(alpha)
// alpha: [0, pi/4]
include <params.scad>
// include <icosphere.scad>
include <egg.scad>

module support_ball_holes() {
  for(i = [0:2])
    rotate([0, 0, 120 * i])
      rotate([-50, 0, 0])
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
  rotate([egg_tilt, 0, 0])
    union() {
      half_track_egg();
      mirror([1, 0, 0])
        half_track_egg();
    }
}

module shell_vtk(path) {
  translate([0, 0, btm_h])
    // rotate([0, 0, 180])
    rotate([0, 180, 0])
      rotate([-90, 0, 0])
        import(path);
}

module shell() {
  // shell_base();
  // shell_3d_minkowski();
  // shell_2d_offset();
  shell_vtk("../../surface-mkw=1.6.stl");
}

btn_offset_x = 30;
btn_offset_y = 33;
btn_ear_thick = 3;

module button_section(extrude_h = 10, offset_r = 0) {
  translate([0, btn_offset_y, 2 - center[2]])
    rotate([8, 0, 0])
      rotate([0, -90, 0])
        rotate([180, 0, 0])
          linear_extrude(extrude_h)
            offset(offset_r)
              scale(18)
                egg_2d(0, 45);
}

switch_d1 = 7.6;
switch_d2 = 2.4;

module switch_hole_section() {
  // switch 6 (body) + 3.4 (pins)
  // pcb = 1.6
  // --> 7.6 + 2.4
  mirror([1, 0, 0])
    union() {
      translate([switch_d1 / 2, 0, 0])
        square([switch_d1, 30.4], center = true);
      translate([switch_d1 + switch_d2 / 2 - 0.01, 0, 0])
        square([switch_d2, 18], center = true);
    }
}

module switch_hole() {
  h = 8.8;
  translate([btn_offset_x - btn_ear_thick + 0.01, btn_offset_y - 3, -center[2] - 0.01]) {
    linear_extrude(h)
      switch_hole_section();
    translate([-switch_d1, +12, h / 2])
      rotate([0, +90, 0])
        cylinder(d = 2.2, h = 8, $fn = 6, center = true);
    translate([-switch_d1, -12, h / 2])
      rotate([0, +90, 0])
        cylinder(d = 2.2, h = 8, $fn = 6, center = true);
  // translate([0, 0, h - 0.01])
  //   linear_extrude(12, scale = 0)
  //     switch_hole_section();
  }
}

module button_hole_left() {
  union() {
    translate([btn_offset_x - 0.01, 0, 0])
      button_section(12, 0.5);
    translate([btn_offset_x, 0, 0]) {
      mirror([1, 0, 0])
        button_section(btn_ear_thick, 2.5);
    }
    switch_hole();
  }
}

module treg_btn() {
  difference() {
    union() {
      hull() {
        btn_r = 1.0;
        for(th = [0:10:90]) {
          dx = btn_r * sin(th);
          translate([dx, 0, 0])
            intersection() {
              shell();
              translate([btn_offset - dx, 0, 0])
                button_section(12, btn_r * (cos(th) - 1));
            }
        }
      }
      translate([btn_offset, 0, 0])
        mirror([1, 0, 0])
          button_section(2.0, 2.0);
    }
    translate([0, 0, -200 - center[2]])
      cube(400, center = true);
  }
}

module treg_top() {
  difference() {
    shell();
    translate([0, 0, -200 - center[2]])
      cube(400, center = true);
    support_ball_holes();
    button_hole_left();
    mirror([1, 0, 0])
      button_hole_left();
  }
}

treg_top();
// translate([13, 0, 0])
//   treg_btn();
// mirror([1, 0, 0])
//   translate([13, 0, 0])
//     treg_btn();

// ball
if (false)
  color("white")
    translate(center)
      scale(ball_r)
        import("icosphere-4.stl");