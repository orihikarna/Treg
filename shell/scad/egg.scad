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

// egg(0, 42, 2);