include <params.scad>

pcba_size = [32, 48];
pcba_screw_pos = 3;
pcba_sensor_pos = 12;

module pcba_hole() {
  translate([0, pcba_size[1] / 2 - pcba_sensor_pos, -center[2] - pcba_offset_z])
    mirror([0, 0, 1])
      linear_extrude(base_h + _clr)
        square(pcba_size, center = true);// pcba plate
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