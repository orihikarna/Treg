include <params.scad>

// switch 6 (body) + 3.4 (pins)
// pcb = 1.6
switch_w = 8.8;
switch_d1 = 8.0;// 6 + 1.6 = 7.6
switch_d2 = 3.0;// 3.4 - 1.6 = 1.8 -> 2.4
switch_h1 = 30.6;
switch_h2 = 18;
switch_angle = 35;
switch_wire_d = switch_d1 + switch_d2;
switch_wire_w = 3;

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
  translate([btn_offset_x - btn_ear_hole_thick, btn_offset_y, -center[2]]) {
    translate([0, 6, 3])
      rotate([switch_angle, 0, 0])
        union() {
          translate([0, 0, -base_h * 0])
            linear_extrude(switch_w + _clr + base_h * 0)
              switch_hole_section();
          for(sgn = [-1, +1])
            translate([-switch_d1, sgn * 12, switch_w / 2])
              rotate([0, +90, 0])
                rotate([0, 0, -switch_angle])
                  cylinder(d = 1.9, h = 8, $fn = 6, center = true);
        }
    // wire hole (verticcal)
    translate([-switch_wire_d / 2, 0, -base_h - _clr])
      linear_extrude(13)
        square([switch_wire_d + _clr, switch_wire_w], center = true);
    // wire hole (horizontal to pcba)
    // w = btn_offset_x - btn_ear_hole_thick - (switch_d1 + switch_d2) - pcba_size[0] / 2;
    w = 5;
    translate([-switch_wire_d - w / 2, 0, -pcba_support_h]) {
      linear_extrude(switch_wire_w / 2 * tan(40), scale = [1, 0])
        square([w + _clr, switch_wire_w], center = true);
      translate([0, 0, _clr])
        mirror([0, 0, 1])
          linear_extrude(base_h + _clr)
            square([w + _clr, switch_wire_w], center = true);
    }
  }
}