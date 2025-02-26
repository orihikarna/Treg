include <icosphere.scad>

// https://math.stackexchange.com/questions/267091/expression-of-the-equations-of-3d-egg-shape-in-terms-of-degrees

// scale(alpha) = 2*cos(alpha) - 1
// integral(1/scale(slapha) * dalpha) = 2/sqrt(3) * tanh^-1(sqrt(3) * tan(alpha / 2))
// https://www.wolframalpha.com/input?i=integration+of+1%2F%282*cos%28x%29-1%29&lang=ja

// tanh(x) = (e^x - e^-x) / (e^x + e^-x) = (e^2x - 1) / (e^2x + 1) = y
// x = artanh(y)
// e^2x = (1 + y) / (1 - y)
function tanh(x) = (exp(2 * x) - 1) / (exp(2 * x) + 1);
function tanh_inv(y) = ln((1 + y) / (1 - y)) / 2;

function L(alpha) = (2 / sqrt(3)) * tanh_inv(sqrt(3) * tan(alpha / 2));
function L_inv(L) = atan(tanh(L * (sqrt(3) / 2)) / sqrt(3)) * 2;

module egg(btm_alpha, top_alpha, recursion = 2) {
  L0 = L(btm_alpha);
  L1 = L(top_alpha);
  hull() {
    for(alpha = [btm_alpha, top_alpha])
      translate([0, 0, tan(alpha)])
        scale((2 - 1 / cos(alpha)))
          icosphere(1, recursion);
    N = round(pow(2, recursion + 2) * 1.4);
    K = round(pow(2, recursion) * 1.5);
    for(k = [1:K - 1]) {
      alpha = L_inv((L1 - L0) * k / K + L0);
      translate([0, 0, 2 * sin(alpha)])
        linear_extrude(0.01, scale = 0)
          rotate([0, 0, 180 / N * k])
            circle(r = 2 * cos(alpha) - 1, $fn = N);
    }
  }
}



// point on ellipse     p = (2*cos(alpha)-1, k*sin(alpha))
function p_x(alpha, k) = 2 * cos(alpha) - 1;
function p_y(alpha, k) = k * 2 * sin(alpha);

// tangential vector    u = (-2*sin(alpha), k*cos(alpha))
// normal vector        n = (-k*cos(alpha), -2*sin(alpha))
function n_x(alpha, k) = -k * 2 * cos(alpha);
function n_y(alpha, k) = -2 * sin(alpha);

// line                 p + n*t = (0, y)
//     p(x) + n(x)*t = 2*cos(alpha) - 1 - k*cos(alpha)*t = 0
function t(alpha, k) = -p_x(alpha, k) / n_x(alpha, k);

// center y = p(y) + n(y)*t
function center_y(alpha, k) = p_y(alpha, k) + n_y(alpha, k) * t(alpha, k);

// radius r = ||normal|| * t
function radius(alpha, k) = sqrt(n_x(alpha, k) * n_x(alpha, k) + n_y(alpha, k) * n_y(alpha, k)) * abs(t(alpha, k));


module egg_eliptic(top_alpha, z_scale = 1, recursion = 2) {
  L = L(top_alpha);
  hull() {
    for(alpha = [0, top_alpha])
      translate([0, 0, center_y(alpha, z_scale)])
        // scale((2 - 1 / cos(alpha)))
        scale(radius(alpha, z_scale))
          icosphere(1, recursion);
    N = round(pow(2, recursion + 2.5));
    K = round(pow(2, recursion) * 1.5 * z_scale);
    for(k = [1:K - 1]) {
      alpha = L_inv(L * k / K);
      translate([0, 0, p_y(alpha, z_scale)])
        linear_extrude(0.01, scale = 0)
          rotate([0, 0, 180 / N * k])
            circle(r = p_x(alpha, z_scale), $fn = N);
    }
  }

}

// translate([-2.2, 0, 0.5])
//   scale([1, 1, 2])
//     egg(0, 43, 3);

// scale([1, 1, 1.5])
//   egg_eliptic(45, 1.5, 3);

// translate([+2.2, 0, 0])
//   egg_eliptic(45, 2.3, 3);