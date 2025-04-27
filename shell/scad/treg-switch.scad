include <params.scad>

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