module pillar(xr, zh) {
  top_h = xr * 1.732;
  translate([0, 0, -0.01])
    linear_extrude(zh - top_h + 0.02)
      circle(r = xr);
  translate([0, 0, zh - top_h])
    linear_extrude(top_h, scale = [0.5 / xr, 1])
      circle(r = xr);
}

module linear_slit(xr, yd, zh1, zh2) {
  hull() {
    translate([0, +yd / 2, -0.01])
      pillar(xr, zh1);
    translate([0, -yd / 2, -0.01])
      pillar(xr, zh2);
  }
}

XR = 7;

difference() {
  import("track-egg.stl");
  if (false) {
    translate([0, 15, 0])
      linear_slit(XR, 40, 45, 47);
    translate([0, -43, 0])
      linear_slit(XR, 40, 42, 19);
    for(xsign = [-1, +1]) {
      translate([xsign * (XR * 2 + 3), 0, 0]) {
        translate([0, +30, 0])
          linear_slit(XR, 20, 39, 43);
        translate([0, -10, 0])
          linear_slit(XR, 20, 44, 42);
        translate([0, -50, 0])
          linear_slit(XR, 20, 34, 16);
      }
      translate([xsign * (XR * 4 + 3 * 2), 0, 0]) {
        translate([0, 24, 0])
          linear_slit(XR, 20, 27, 31);
        translate([0, -16, 0])
          linear_slit(XR, 20, 31, 25);
      }
    }
  }
  if (true) {
    translate([0, -6, -0.01])
      linear_slit(30, 36, 42, 34);
  }
//   translate([-200, 0, 0])
//   translate([-200 + (XR * 2 + 3), 0, 0])
//   translate([-200 + (XR * 4 + 3 * 2), 0, 0])
// cube(400, center = true);
}