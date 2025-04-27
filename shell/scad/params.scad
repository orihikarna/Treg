_clr = 0.01;

ball_r = 57.2 / 2;
hole_r = ball_r + 2.4 / 2;
btm_h = 0 + 0;
base_h = 9;// + 1;

center = [0, -30, hole_r + btm_h];
egg_btm_alpha = 0;
egg_top_alpha = 36;
egg_scale = [38, 42, 68];
egg_tilt = 12;

mkw_r = 1.6;

hole_mkw_r = hole_r + mkw_r;
egg_mkw_scale = [egg_scale[0] - mkw_r, egg_scale[1] - mkw_r, egg_scale[2] - mkw_r];

btn_offset_x = 31;
btn_offset_y = 24;
btn_ear_hole_gap = 0.3;
btn_ear_thick = 1.6;
btn_ear_btm_thick = 1.2;
btn_ear_hole_thick = 2.6;
btn_ear_roffset = 1.6;
btn_ear_hole_roffset = btn_ear_roffset + btn_ear_hole_gap + 0.1;