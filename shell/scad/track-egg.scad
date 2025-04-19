include <track-shell.scad>


btn_offset_x = 31;
btn_offset_y = 24;
btn_ear_hole_gap = 0.3;
btn_ear_thick = 1.6;
btn_ear_hole_thick = 2.6;
btn_ear_roffset = 1.6;
btn_ear_hole_roffset = btn_ear_roffset + btn_ear_hole_gap + 0.1;

switch_d1 = 8.0;
switch_d2 = 3.0;
switch_w = 8.8;
switch_h1 = 30.6;
switch_h2 = 18;

module switch_hole_section() {
  // switch 6 (body) + 3.4 (pins)
  // pcb = 1.6
  // --> 7.6 + 2.4
  mirror([1, 0, 0])
    union() {
      translate([switch_d1 / 2 - _clr, 0, 0])
        square([switch_d1 + _clr * 2, switch_h1], center = true);
      translate([switch_d1 + switch_d2 / 2 - _clr, 0, 0])
        square([switch_d2 + _clr * 2, switch_h2], center = true);
    }
}

module switch_hole() {
  // translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y + 2, -center[2] - _clr]) {
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y + 5, 2 - center[2] - _clr]) {
    rotate([45, 0, 0])
      union() {
        translate([0, 0, -base_h * 0])
          linear_extrude(switch_w + _clr + base_h * 0)
            switch_hole_section();
        // translate([0, 0, switch_w])
        //   linear_extrude(12, scale = 0.1)
        //     switch_hole_section();
        for(sgn = [-1, +1])
          translate([-switch_d1, sgn * 12, switch_w / 2])
            rotate([0, +90, 0])
              cylinder(d = 1.9, h = 8, $fn = 6, center = true);
      }
  }
}

module button_section(extrude_h = 10, offset_r = 0) {
  translate([0, btn_offset_y, -2 - center[2]])
    rotate([60, 0, 0])
      rotate([0, +90, 0])
        linear_extrude(extrude_h)
          offset(offset_r)
            scale(20)
              egg_2d(0, 45);
}

spt_w = 3;
spt_dy = 15;

module button_support_hole() {
  x1 = spt_w - 0.4;
  x2 = spt_w * 2 - x1;
  y = spt_w + 0.8;
  h1 = spt_w + 1;
  h2 = spt_w * 2 + 2;
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y, -center[2] - base_h])
    for(sgn_y = [-1, +1])
      translate([0, spt_dy * sgn_y, 0])
        mirror([1, 0, 0]) {
          translate([x1 / 2, 0, -_clr]) {
            linear_extrude(h1 + _clr)
              square([x1 + _clr, y], center = true);
            translate([0, 0, h1])
              linear_extrude(y * 0.4, scale = [1, 0])
                square([x1 + _clr, y], center = true);
          }
          translate([x1 + x2 / 2, 0, -_clr]) {
            linear_extrude(h2 + _clr)
              square([x2 + _clr, y], center = true);
            translate([0, 0, h2])
              linear_extrude(y * 0.4, scale = [1, 0])
                square([x2 + _clr, y], center = true);
          }
        }
}

module button_support() {
  x1 = spt_w;
  x2 = spt_w * 2 - x1;
  y = spt_w;
  h1 = spt_w;
  h2 = spt_w * 2;
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y, -center[2] - base_h])
    for(sgn_y = [-1, +1])
      translate([0, spt_dy * sgn_y, 0])
        mirror([1, 0, 0]) {
          translate([(1.0 + x1 + x2 / 2) / 2 - 0.5, 0, -_clr]) {
            linear_extrude(h1 + _clr)
              square([1.0 + x1 + x2 / 2 + _clr, y], center = true);
          }
          translate([x1 + x2 / 2, 0, -_clr]) {
            linear_extrude(h2 + _clr, scale = [0.7, 1])
              square([x2 + _clr, y], center = true);
          }
        }
}

module button_ear_hole() {
  union() {
    // outer side
    translate([btn_offset_x - _clr, 0, 0])
      button_section(16, btn_ear_hole_gap);
    // innder side
    translate([btn_offset_x, 0, 0]) {
      mirror([1, 0, 0])
        hull() {
          d = btn_ear_roffset * 0.8;// same as in treg_btn()
          button_section(0.1, btn_ear_hole_gap);
          translate([d, 0, 0])
            button_section(btn_ear_hole_thick - d, btn_ear_hole_roffset);
        }
    }
  }
}

module button_hole() {
  union() {
    switch_hole();
    button_ear_hole();
    button_support_hole();
    if (true)
      translate([0, 0, -center[2] - base_h - _clr])
        union() {
          linear_extrude(base_h - 2 + _clr)
            intersection() {
              translate([btn_offset_x, btn_offset_y - 20, 0])
                square(40, center = true);
              projection(cut = false)
                button_ear_hole();
            }
          linear_extrude(base_h + 7 + _clr)
            intersection() {
              translate([btn_offset_x, btn_offset_y + 20, 0])
                square(40, center = true);
              projection(cut = false)
                button_ear_hole();
            }
        }
  }
}

module treg_btn_base() {
  difference() {
    union() {
      hull() {
        btn_r = 1.6;
        for(th = [0:10:90]) {
          dx = btn_r * sin(th);
          translate([dx, 0, 0])
            intersection() {
              shell();
              translate([btn_offset_x - dx, 0, 0])
                button_section(16, btn_r * (cos(th) - 1));
            }
        }
      }
      d = btn_ear_roffset * 0.8;
      translate([btn_offset_x, 0, 0])
        mirror([1, 0, 0]) {
          hull() {
            button_section(0.1, 0);
            translate([d, 0, 0])
              button_section(btn_ear_thick - d, btn_ear_roffset);
          }
        }
    }
    translate([0, 0, -200 - center[2] - base_h])
      cube(400, center = true);
  }
}

module treg_btn() {
  hinge_offset_x = -0.2;
  hinge_offset_y = 0.25;
  hinge_offset_z = 2.4;
  hinge_outer_r = 2.5;
  hinge_size_y = 43.5;
  difference() {
    union() {
      treg_btn_base();
      translate([btn_offset_x - btn_ear_thick + hinge_offset_x, btn_offset_y, -center[2] - base_h + hinge_offset_z]) {
        // button_support();
        difference() {
          translate([0, hinge_offset_y, 0])
            union() {
              rotate([90, 0, 0])
                cylinder(r = hinge_outer_r, h = hinge_size_y, center = true, $fn = 60);
              // translate([0, 0, -5])
              // cube([hinge_outer_r * 2, hinge_size_y, 10], center = true);
              translate([-0, 0, -hinge_outer_r * 0.5])
                mirror([0, 0, 1])
                  linear_extrude(hinge_outer_r * tan(60), scale = [0, 1])
                    square([hinge_outer_r * 2 * cos(30), hinge_size_y], center = true);
            }
          cube([10, 30, 10], center = true);
          translate([25 - hinge_offset_x + _clr, 0, 0])
            cube(50, center = true);
        }
      }
    }
    translate([btn_offset_x - btn_ear_thick + hinge_offset_x, btn_offset_y, -center[2] - base_h + hinge_offset_z]) {
      // button_support();
      rotate([0, 30, 0])
        rotate([90, 0, 0])
          cylinder(d = 2.2, h = 48, center = true, $fn = 6);
      translate([0, 0, -5])
        cube([2.2 * cos(30), 30, 10], center = true);
    }
  }
}

pcba_size = [32, 48];
pcba_screw_pos = 3;
pcba_sensor_pos = 12;

module pcba_hole() {
  translate([0, pcba_size[1] / 2 - pcba_sensor_pos, -center[2] - 0.6])
    mirror([0, 0, 1]) {
      linear_extrude(base_h + _clr) {
        square(pcba_size, center = true);
        translate([0, 12.4, 0])
          square([100, 3], center = true);
      }
    }
  translate([0, 0, -center[2]])
    cylinder(d = 8, h = 2, center = true, $fn = 90);
}

module pcba_screw_support() {
  translate([0, pcba_size[1] / 2 - pcba_sensor_pos, -center[2]])
    mirror([0, 0, 1]) {
      for(y = [-1, +1])
        for(x = [-1, +1]) {
          translate([(pcba_size[0] / 2 - pcba_screw_pos) * x, (pcba_size[1] / 2 - pcba_screw_pos) * y, 0])
            difference() {
              r = 3;
              h = 7.4 - 1.2 - 1.6;
              translate([0, 0, h / 2])
                union() {
                  cylinder(h = h, r = r, center = true, $fn = 90);
                  translate([r / 2 * x, 0, 0])
                    cube([r, 2 * r, h], center = true);
                  translate([0, r / 2 * y, 0])
                    cube([2 * r, r, h], center = true);
                }
              translate([0, 0, h - 4 / 2 - 0.2])
                cylinder(h = 4, d = 2.0, center = true, $fn = 6);
            }
        }
    }
}

module treg_top() {
  difference() {
    intersection() {
      shell();
      translate([0, -200 + 60, 0])
        cube(400, center = true);
    }
    support_ball_holes();
    translate([0, 0, -200 - center[2] - base_h])
      cube(400, center = true);
    button_hole();
    mirror([1, 0, 0])
      button_hole();
    pcba_hole();
  }
  pcba_screw_support();
}

if (false) {
  intersection() {
    treg_top();
  // translate([0, 0, -200 - 20])
  //   cube(400, center = true);
  }
} else {
  intersection() {
    // rotate([0, 2, 0])
    translate([-btn_offset_x + btn_ear_thick + spt_w * 2, 0, center[2] + base_h])
      treg_btn();
    translate([0, 0, 200])
      cube(400, center = true);
  }
}

// ball
if (false)
  color("white")
    translate(center)
      scale(ball_r)
        import("icosphere-4.stl");