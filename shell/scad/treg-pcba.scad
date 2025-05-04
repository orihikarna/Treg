include <params.scad>

pcba_screw_pos = 3;
pcba_sensor_pos = 12;

module pcba_hole() {
  border_width = pcba_support_h * 2 * tan(90 - 40);
  difference() {
    translate([0, pcba_bare_size[1] / 2 - pcba_sensor_pos, -center[2] - pcba_offset_z])
      union() {
        mirror([0, 0, 1])
          linear_extrude(base_h + _clr)
            square(pcba_hole_size, center = true);// pcba plate
        roof_x = pcba_hole_size[0] - border_width;
        roof_y = (pcba_hole_size[1] - border_width) * 0.4;
        translate([0, (pcba_hole_size[1] - border_width) / 2 - roof_y / 2, -_clr])
          linear_extrude(roof_y / 2 * tan(40), scale = [0.5, 0])
            square([roof_x, roof_y], center = true);
      }
    translate([0, pcba_sensor_pos, -center[2] - pcba_offset_z])// boarder slopes
      mirror([0, 0, 1]) {
        for(x = [-1, +1])
          translate([pcba_hole_size[0] / 2 * x, 0, 0])
            linear_extrude(pcba_support_h, scale = [0, 1])
              square([border_width, pcba_hole_size[1]], center = true);
        for(y = [-1, +1])
          translate([0, pcba_hole_size[1] / 2 * y, 0])
            linear_extrude(pcba_support_h, scale = [1, 0])
              square([pcba_hole_size[0], border_width], center = true);
      }
  }
  translate([0, 0, -center[2]])
    cylinder(d = 8, h = 2, center = true, $fn = 90);// ball sensor
}

module pcba_screw_support() {
  h_add = 3;
  translate([0, pcba_bare_size[1] / 2 - pcba_sensor_pos, -center[2] - pcba_offset_z])
    mirror([0, 0, 1]) {
      for(y = [-1, +1])
        for(x = [-1, +1]) {
          translate([(pcba_bare_size[0] / 2 - pcba_screw_pos) * x, (pcba_bare_size[1] / 2 - pcba_screw_pos) * y, 0])
            difference() {
              r = 3.3;
              translate([0, 0, pcba_support_h / 2])
                union() {
                  translate([0, 0, -h_add / 2])
                    cylinder(h = pcba_support_h + h_add, r = r, center = true, $fn = 90);
                  translate([r / 2 * x, 0, 0])
                    cube([r, 2 * r, pcba_support_h], center = true);
                  translate([0, r / 2 * y, 0])
                    cube([2 * r, r, pcba_support_h], center = true);
                }
              translate([0, 0, pcba_support_h - 4 / 2 - 0.2])
                cylinder(h = 4, d = 2.0, center = true, $fn = 6);
            }
        }
    }
}