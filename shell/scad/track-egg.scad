include <track-shell.scad>
include <treg-switch.scad>
include <treg-button.scad>
include <treg-pcba.scad>

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