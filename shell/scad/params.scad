ball_r = 57.2 / 2;
hole_r = ball_r + 2.5 / 2;
btm_h = 1;
center = [0, -36, hole_r + btm_h];
egg_btm_alpha = 0;
egg_top_alpha = 36;
egg_scale = [35, 40, 80];

mkw_r = 3;

hole_mkw_r = hole_r + mkw_r;
egg_mkw_scale = [egg_scale[0] - mkw_r, egg_scale[1] - mkw_r, egg_scale[2] - mkw_r];