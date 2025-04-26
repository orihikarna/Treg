include <track-shell.scad>

// switch 6 (body) + 3.4 (pins)
// pcb = 1.6
switch_w = 8.8;
switch_d1 = 8.0;// 6 + 1.6 = 7.6
switch_d2 = 3.0;// 3.4 - 1.6 = 1.8 -> 2.4
switch_h1 = 30.6;
switch_h2 = 18;

module switch_hole_section() {
  mirror([1, 0, 0])
    union() {
      translate([switch_d1 / 2 - _clr, 0, 0])
        square([switch_d1 + _clr * 2, switch_h1], center = true);
      translate([switch_d1 + switch_d2 / 2 - _clr, 0, 0])
        square([switch_d2 + _clr * 2, switch_h2], center = true);
    }
}

module switch_hole() {
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y + 5, 2 - center[2] - _clr]) {
    rotate([45, 0, 0])
      union() {
        translate([0, 0, -base_h * 0])
          linear_extrude(switch_w + _clr + base_h * 0)
            switch_hole_section();
        for(sgn = [-1, +1])
          translate([-switch_d1, sgn * 12, switch_w / 2])
            rotate([0, +90, 0])
              cylinder(d = 1.9, h = 8, $fn = 6, center = true);
      }
  }
}

btn_offset_x = 31;
btn_offset_y = 24;
btn_ear_hole_gap = 0.3;
btn_ear_thick = 1.6;
btn_ear_btm_thick = 1.2;
btn_ear_hole_thick = 2.6;
btn_ear_roffset = 1.6;
btn_ear_hole_roffset = btn_ear_roffset + btn_ear_hole_gap + 0.1;

module button_section(extrude_h = 10, offset_r = 0) {
  translate([0, btn_offset_y, -2 - center[2]])
    rotate([60, 0, 0])
      rotate([0, +90, 0])
        linear_extrude(extrude_h)
          offset(offset_r)
            scale(20)
              egg_2d(0, 45);
}

module button_shape() {
  hull()
    for(th = [0:10:90]) {
      btn_r = 1.6;
      dx = btn_r * sin(th);
      offset_r = btn_r * (cos(th) - 1);
      translate([dx, 0, 0])
        intersection() {
          shell();
          translate([btn_offset_x - dx, 0, 0])
            button_section(16, offset_r);
        }
    }
}

module button_hole() {
  union() {
    button_ear_hole();
    // button_hook_support_hole();
    button_hinge_support_hole();
    if (true)// button hole extruded downto bottom plate
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

module button_ear() {
  translate([btn_offset_x, 0, 0])
    mirror([1, 0, 0])
      hull() {
        button_section(0.1, 0);
        translate([btn_ear_btm_thick, 0, 0])
          button_section(btn_ear_thick - btn_ear_btm_thick, btn_ear_roffset);
      }
}

module button_ear_hole() {
  union() {
    // outer side
    translate([btn_offset_x - _clr, 0, 0])
      button_section(16, btn_ear_hole_gap);
    // inner side
    translate([btn_offset_x, 0, 0])
      mirror([1, 0, 0])
        hull() {
          button_section(0.1, btn_ear_hole_gap);
          translate([btn_ear_btm_thick, 0, 0])
            button_section(btn_ear_hole_thick - btn_ear_btm_thick, btn_ear_hole_roffset);
        }
  }
}

hook_w = 3;
hook_dy = 15;

module button_hook_support_hole() {
  x1 = hook_w - 0.4;
  x2 = hook_w * 2 - x1;
  y = hook_w + 0.8;
  h1 = hook_w + 1;
  h2 = hook_w * 2 + 2;
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y, -center[2] - base_h])
    for(sgn_y = [-1, +1])
      translate([0, hook_dy * sgn_y, 0])
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

module button_hook_support() {
  x1 = hook_w;
  x2 = hook_w * 2 - x1;
  y = hook_w;
  h1 = hook_w;
  h2 = hook_w * 2;
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y, -center[2] - base_h])
    for(sgn_y = [-1, +1])
      translate([0, hook_dy * sgn_y, 0])
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

hinge_offset_x = -0.2;
hinge_offset_y = 0.25;
hinge_offset_z = 2.4;
hinge_outer_r = 2.5;
hinge_size_y = 43.5;
hinge_size_y_inner = 30;

module button_hinge_support() {
  angle = 10;
  translate([btn_offset_x - btn_ear_thick + hinge_offset_x, btn_offset_y, -center[2] - base_h + hinge_offset_z])
    difference() {
      translate([0, hinge_offset_y, 0])
        union() {
          rotate([90, 0, 0])
            cylinder(r = hinge_outer_r, h = hinge_size_y, center = true, $fn = 60);
          translate([0, 0, -hinge_outer_r * sin(angle)])
            mirror([0, 0, 1])
              linear_extrude(hinge_outer_r * cos(angle) * tan(90 - angle), scale = [0, 1])
                square([hinge_outer_r * 2 * cos(angle), hinge_size_y], center = true);
        }
      cube([10, hinge_size_y_inner, 10], center = true);
      translate([25 - hinge_offset_x + _clr, 0, 0])
        cube(50, center = true);
    }
}

d_pin_hole = 2.00;

module button_hinge_support_hole_button() {
  translate([btn_offset_x - btn_ear_thick + hinge_offset_x, btn_offset_y, -center[2] - base_h + hinge_offset_z]) {
    rotate([0, 30, 0])
      rotate([90, 0, 0])
        cylinder(d = d_pin_hole / cos(30), h = 48, center = true, $fn = 6);
    translate([0, 0, -5])
      cube([2.2 * cos(30), 30, 10], center = true);
  }
}

module button_hinge_support_hole() {
  gap = 0.5;
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y, -center[2] - base_h + hinge_offset_z])
    difference() {
      translate([0, hinge_offset_y, 0])
        union() {
          rotate([90, 0, 0])
            cylinder(r = hinge_outer_r + gap, h = hinge_size_y + gap * 2, center = true, $fn = 60);
          translate([0, 0, (hinge_outer_r + gap) * sin(50)])
            linear_extrude((hinge_outer_r + gap) * cos(50) * tan(40), scale = [0, 1])
              square([(hinge_outer_r + gap) * 2 * cos(50), hinge_size_y + gap * 2], center = true);
          mirror([0, 0, 1])
            linear_extrude(hinge_outer_r * 1.2)
              square([(hinge_outer_r + gap) * 2, hinge_size_y + gap * 2], center = true);
          translate([0, 1.5, 0])
            rotate([0, 30, 0])
              rotate([90, 0, 0])
                cylinder(d = 1.85 * 2 / sqrt(3), h = 60, center = true, $fn = 6);
        }
      difference() {
        cube([10, hinge_size_y_inner - gap * 2, 10], center = true);
        rotate([0, 30, 0])
          rotate([90, 0, 0])
            cylinder(d = d_pin_hole, h = 48, center = true, $fn = 6);
        translate([0, 0, -5])
          cube([d_pin_hole * cos(30), hinge_size_y_inner - gap * 2, 10], center = true);
      }
      translate([5 + d_pin_hole / 2 * cos(30), 0, 0])
        cube([10, 50, 10], center = true);
    }
}

module treg_button() {
  difference() {
    union() {
      button_shape();
      button_ear();
      // button_hook_support();
      button_hinge_support();
    }
    button_hinge_support_hole_button();
    translate([0, 0, -200 - center[2] - base_h])
      cube(400, center = true);
  }
}

pcba_size = [32, 48];
pcba_screw_pos = 3;
pcba_sensor_pos = 12;

module pcba_hole() {
  translate([0, pcba_size[1] / 2 - pcba_sensor_pos, -center[2] - 0.6])
    mirror([0, 0, 1]) {
      linear_extrude(base_h + _clr)
        square(pcba_size, center = true);// pcba plate
      translate([0, 12.4, 3]) {// switch wire
        linear_extrude(base_h)
          square([100, 3], center = true);
        translate([0, 0, _clr])
          mirror([0, 0, 1])
            linear_extrude(3.0 / 2, scale = [1, 0]) {
              square([100, 3], center = true);
            }
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
    shell();
    support_ball_holes();
    switch_hole();
    button_hole();
    mirror([1, 0, 0]) {
      switch_hole();
      button_hole();
    }
    pcba_hole();
    translate([0, 0, -200 - center[2] - base_h])
      cube(400, center = true);// cut bottom plate
  }
  pcba_screw_support();
}

if (true) {
  intersection() {
    treg_top();
    // translate([0, 0, -200 - 20])
    //   cube(400, center = true);
    translate([0, -200 + 60, 0])
      cube(400, center = true);// extract ball & button part
    translate([btn_offset_x, btn_offset_y, -center[2] - base_h + 22])
      cube([30, 70, 44], center = true);// extract button area
  }
}
if (false) {
  intersection() {
    // rotate([0, 2, 0])
    translate([btn_ear_thick - btn_ear_hole_thick, 0, 0])
      treg_button();
    translate([0, 0, -center[2] - base_h + 0.2 + 200])
      cube(400, center = true);
  }
}
if (false)
  color("blue", 0.8)
    button_hinge_support_hole();