#include "TrackEgg.h"

namespace {
const Eigen::Vector3f ball_ctr{0, ball_y, ball_z};
const Eigen::Vector3f EGG_OOB = Eigen::Vector3f::Zero();
const TrackEggRidge ridge(720);
}  // namespace

Eigen::Vector3f egg_surface(float z, float theta) {
  float ez = (z - egg_org_z) / egg_scale_z;
  float r = 0;
  if (ez < egg_zmin) {
    return EGG_OOB;
  } else if (ez < egg_zedge_btm) {
    const float dz = ez - egg_zbtm;
    r = std::sqrt(egg_rad_btm * egg_rad_btm - dz * dz);
  } else if (ez < egg_zedge_top) {
    const float _alpha = std::asin(ez / 2);
    r = 2 * std::cos(_alpha) - 1;
  } else if (ez < egg_zmax) {
    const float dz = ez - egg_ztop;
    r = std::sqrt(egg_rad_top * egg_rad_top - dz * dz);
  } else {
    return EGG_OOB;
  }
  const float zz = ez * egg_scale_z + egg_org_z;
  const float xx = std::cos(theta) * r * egg_scale_x;
  const float yy = std::sin(theta) * r * egg_scale_y;
  return {xx, yy, zz};
}

Eigen::Vector3f hole_surface(float azim, float elev, float r = hole_r) {
  Eigen::Vector3f pos;
  pos[2] = (r * std::cos(elev)) * std::cos(azim);
  pos[0] = (r * std::cos(elev)) * std::sin(azim);
  pos[1] = (r * std::sin(elev));
  return pos + ball_ctr;
}

std::tuple<float, float> hole_pos2azel(const Eigen::Vector3f &_pos) {
  const Eigen::Vector3f pos = _pos - ball_ctr;
  const float r = std::sqrt(pos[2] * pos[2] + pos[0] * pos[0]);
  const float elev = std::atan2(pos[1], r);
  const float azim = std::atan2(pos[0], pos[2]);
  return {azim, elev};
}

bool is_inside_hole(const Eigen::Vector3f &pos) {
  constexpr float hole_r2 = hole_r * hole_r;
  return (pos - ball_ctr).squaredNorm() <= hole_r2;
}

// ===== 2. distance to egg surface for Fast Marching =====
class DistanceToEggSurface : public vnl_cost_function {
 private:
  Eigen::Vector3f pos_;

 public:
  DistanceToEggSurface(const Eigen::Vector3f &pos) : pos_(pos), vnl_cost_function(2) {}

  double f(vnl_vector<double> const &x) override {
    const float z = x[0];
    const float theta = x[1];
    const Eigen::Vector3f pos_egg = egg_surface(z, theta);
    if (pos_egg == EGG_OOB) return 1e9;
    const float dist2 = (pos_ - pos_egg).squaredNorm();
    return dist2;
  }
};

std::tuple<float, Eigen::Vector3f> calc_dist_egg(const Eigen::Vector3f &pos) {
  DistanceToEggSurface dist(pos);
  vnl_powell minimizer(&dist);
  minimizer.set_f_tolerance(1e-6);
  minimizer.set_trace(true);
  double z = pos[2];
  {  // bring it inside the egg
    z = std::max(z, egg_zmin + egg_org_z + 1.0);
    z = std::min(z, egg_zmax + egg_org_z - 1.0);
  }
  vnl_vector<double> x(2);
  x[0] = z;
  x[1] = atan2(pos[1], pos[0]);
  minimizer.minimize(x);
  const Eigen::Vector3f pos_egg = egg_surface(x[0], x[1]);
  // distance sign
  const bool inside = ((pos - pos_egg).dot(pos_egg - Eigen::Vector3f(0, 0, egg_org_z)) < 0);
  const float dist2 = minimizer.get_end_error();
  const float _dist = std::sqrt(dist2) * ((inside) ? -1 : +1);
  return {_dist, pos_egg};
}

std::tuple<float, Eigen::Vector3f> calc_dist_hole(const Eigen::Vector3f &_pos) {
  const Eigen::Vector3f pos = _pos - ball_ctr;
  const float pos_norm = pos.norm();
  const float dist = hole_r - pos_norm;  // inside hole --> positive, outside hole --> negative
  const Eigen::Vector3f pos_hole = pos * (hole_r / pos_norm) + ball_ctr;
  return {dist, pos_hole};
}

class DistanceToEggSurface_BallAzim_1pass : public vnl_cost_function {
 private:
  float azim_;

 public:
  DistanceToEggSurface_BallAzim_1pass(float azim) : azim_(azim), vnl_cost_function(3) {}

  double f(vnl_vector<double> const &x) override {
    const float elev = x[0];
    const float z = x[1];
    const float theta = x[2];
    if (elev < -M_PI / 2 || +M_PI / 2 < elev) return 1e9;
    const Eigen::Vector3f pos_egg = egg_surface(z, theta);
    if (pos_egg == EGG_OOB) return 1e9;
    const Eigen::Vector3f pos_ball = hole_surface(azim_, elev);
    const float dist2 = (pos_ball - pos_egg).squaredNorm();
    return dist2;
  }
};

float calc_ball_elev_1pass(float azim) {
  DistanceToEggSurface_BallAzim_1pass dist(azim);
  vnl_powell minimizer(&dist);
  minimizer.set_f_tolerance(1e-8);
  minimizer.set_trace(true);
  const float elev0 = -15 * (M_PI / 180);
  const Eigen::Vector3f pos_ball0 = hole_surface(azim, elev0);
  vnl_vector<double> x(3);
  x[0] = elev0;
  x[1] = pos_ball0[2];  // z
  x[2] = 0;             // theta
  minimizer.minimize(x);
  const float dist2 = minimizer.get_end_error();
  if (dist2 > 1e-5) std::cerr << "[ERROR] dist2 = " << dist2 << std::endl;
  return x[0];
}

class DistanceToEggSurface_BallAzim_2pass : public vnl_cost_function {
 private:
  float azim_;

 public:
  DistanceToEggSurface_BallAzim_2pass(float azim) : azim_(azim), vnl_cost_function(1) {}

  double f(vnl_vector<double> const &x) override {
    const float elev = x[0];
    if (elev < -M_PI / 2 || +M_PI / 2 < elev) return 1e9;
    const Eigen::Vector3f pos_ball = hole_surface(azim_, elev);
    Eigen::Vector3f pos_egg;
    {
      DistanceToEggSurface dist(pos_ball);
      vnl_powell minimizer(&dist);
      minimizer.set_f_tolerance(1e-6);
      minimizer.set_trace(true);
      vnl_vector<double> _x(2);
      _x[0] = pos_ball[2];
      _x[1] = atan2(pos_ball[1], pos_ball[0]);
      minimizer.minimize(_x);
      pos_egg = egg_surface(_x[0], _x[1]);
    }
    const float dist2 = (pos_ball - pos_egg).squaredNorm();
    return dist2;
  }
};

float calc_ball_elev_2pass(float azim) {
  DistanceToEggSurface_BallAzim_2pass dist(azim);
  vnl_powell minimizer(&dist);
  minimizer.set_f_tolerance(1e-8);
  minimizer.set_trace(true);
  const float elev0 = -15 * (M_PI / 180);
  vnl_vector<double> x(1);
  x[0] = elev0;
  minimizer.minimize(x);
  const float dist2 = minimizer.get_end_error();
  if (dist2 > 1e-5) std::cerr << "[ERROR] dist2 = " << dist2 << std::endl;
  return x[0];
}

void TrackEggRidge::Init() {
  azims_.resize(N);
  elevs_.resize(N);
  points_.resize(N);
  for (size_t n = 0; n < N; ++n) {
    const int m = n - M;
    const T azim = dazim * m;
    // const T elev = calc_ball_elev_1pass(azim);
    const T elev = calc_ball_elev_2pass(azim);
    azims_[n] = azim;
    elevs_[n] = elev;
    points_[n] = hole_surface(azim, elev);
    // std::cout << "(azim, elev) = " << azim * (180 / M_PI) << ", " << elev * (180 / M_PI) << std::endl;
  }
}

// on the hole surface?
bool TrackEggRidge::IsOnHoleSurface(float azim, float elev) const {
  if (azim < -M_PI || +M_PI < azim) {
    std::cerr << "(" << __LINE__ << ") azim = " << azim << " is out of range [-M_PI, +M_PI]" << std::endl;
    throw;
  }
  const int n0 = int(std::floor(azim / dazim)) + M;
  float elev0 = -101;
  for (int n = std::max(n0 - 1, 0); n < std::min<int>(n0 + 1, N); ++n) {
    if (azims_[n] <= azim && azim <= azims_[n + 1]) {
      const float k = (azim - azims_[n0]) / (azims_[n0 + 1] - azims_[n0]);
      elev0 = elevs_[n0] + (elevs_[n0 + 1] - elevs_[n0]) * k;
      break;
    }
  }
  if (elev0 < -100) {
    std::cerr << "(" << __LINE__ << ") azim  = " << azim << " not handled" << std::endl;
    throw;
  }
  return (elev < elev0);
}

TrackEggRidge::T TrackEggRidge::CalcMinDist(const Eigen::Vector3f &pos) const {
  T min_dist2 = std::numeric_limits<T>::max();
  for (size_t n = 0; n < N; ++n) {
    const T dist2 = (pos - points_[n]).squaredNorm();
    if (min_dist2 > dist2) {
      min_dist2 = dist2;
    }
  }
  return std::sqrt(min_dist2);
}

std::tuple<NodeContainer::Pointer, NodeContainer::Pointer> TrackEggSeeds() {
  const float min_scale = std::min(std::min(egg_scale_x, egg_scale_y), egg_scale_z);
  const float max_scale = std::max(std::max(egg_scale_x, egg_scale_y), egg_scale_z);

  auto seeds = NodeContainer::New();
  auto outside = NodeContainer::New();
  seeds->Initialize();
  outside->Initialize();
  size_t cnt_seeds = 0;
  size_t cnt_outside = 0;
  InternalImageType::IndexType pos;
  for (int iz = 0; iz < SizeZ; ++iz) {
    pos[2] = iz;
    const float z = spacing * (iz - SizeZ / 2);
    const float sz = z * (1 / egg_scale_z);         // scaled for egg
    const float ez = sz - egg_org_z / egg_scale_z;  // egg
    const float ez2 = ez * ez;                      // egg^2
    const float bz = z - ball_z;                    // ball
    const float bz2 = bz * bz;                      // ball^2
    for (int iy = 0; iy < SizeY; ++iy) {
      pos[1] = iy;
      const float y = spacing * (iy - SizeY / 2);
      const float ey = y * (1 / egg_scale_y);  // scaled for egg
      const float ey2 = ey * ey;               // egg^2
      const float by = y - ball_y;             // ball
      const float by2 = by * by;               // ball^2
      for (int ix = 0; ix < SizeX; ++ix) {
        pos[0] = ix;
        const float x = spacing * (ix - SizeX / 2);
        const float ex = x * (1 / egg_scale_x);  // scaled for egg
        const float ex2 = ex * ex;               // egg^2
        const float bx = x;                      // ball
        const float bx2 = bx * bx;               // ball^2

        // egg
        const float r2 = ex2 + ey2;
        float approx_dist_egg = std::numeric_limits<float>::max();
        if (ez <= 0) {
          approx_dist_egg = std::sqrt(r2 + ez2) - 1;
        } else {
          const float R = std::sqrt(r2) + 1;
          if (std::atan2(ez, R) < egg_alpha) {
            approx_dist_egg = std::sqrt(R * R + ez2) - 2;
          } else {
            const float dz_top = ez - egg_ztop;
            approx_dist_egg = std::sqrt(r2 + dz_top * dz_top) - egg_rad_top;
          }
        }
        const float approx_dist_egg_pxl = approx_dist_egg * min_scale / spacing;  // possible min distance
        // inside egg negative, outside egg positive

        // ball
        const float dist_ball = hole_r - std::sqrt(bx2 + by2 + bz2);  // inside positive, outside negative
        const float dist_ball_pxl = dist_ball / spacing;

        /*
        0. 稜線を求める
        1. egg_surface への最短位置を求める hole の外なら有効（最短距離＝ egg 面に垂直）
        2. hole_surface への最短位置を求める。稜線より下なら有効（最短距離＝ ball 面に垂直）
        3-a. 1か2が両方も有効なら短い方を採用
        3-b. 1か2の片方が有効ならそれを採用？
        3-c, いずれも無効なら、稜線への最短位置を求める

        片方のみが有効の時に、稜線を確認する必要があるか？
        egg 面の外側は凸面なのでない
        hole 面の外側も凸面なのでない
        */
        if (approx_dist_egg_pxl < 2 && dist_ball_pxl < 2) {
          NodeType node;
          node.SetIndex(pos);
          if (approx_dist_egg_pxl < -3.6 && dist_ball_pxl < -3.6) {
            // do nothing
          } else if (approx_dist_egg_pxl < -1.8 && dist_ball_pxl < -1.8) {
            // outside
            outside->InsertElement(cnt_outside++, node);
          } else {
            const Eigen::Vector3f pos{x, y, z};

            const auto [dist_egg, egg_foot] = calc_dist_egg(pos);
            const float dist_egg_pxl = dist_egg / spacing;
            const bool egg_foot_valid = (is_inside_hole(egg_foot) == false);

            const auto [dist_hole, hole_foot] = calc_dist_hole(pos);
            const float dist_hole_pxl = dist_hole / spacing;
            const auto [azim, elev] = hole_pos2azel(pos);
            const bool hole_foot_valid = ridge.IsOnHoleSurface(azim, elev);

            float dist = std::numeric_limits<float>::max();
            if (dist_egg_pxl < 0 && dist_hole_pxl < 0) {  // inside track-egg
              dist = -std::numeric_limits<float>::max();
              if (egg_foot_valid) dist = std::max(dist, dist_egg_pxl);
              if (hole_foot_valid) dist = std::max(dist, dist_hole_pxl);
              if (egg_foot_valid == false) {
                // std::cerr << "(" << __LINE__ << ") egg_foot_valid = false" << std::endl;
                // throw;
              }
              if (hole_foot_valid == false) {
                // std::cerr << "(" << __LINE__ << ") hole_foot_valid = false" << std::endl;
                // throw;
              }
            } else {  // outside track_egg
              if (dist_egg_pxl >= 0 && egg_foot_valid) dist = std::min(dist, dist_egg_pxl);
              if (dist_hole_pxl >= 0 && hole_foot_valid) dist = std::min(dist, dist_hole_pxl);
              dist = std::min(dist, ridge.CalcMinDist(pos) / spacing);
            }
            node.SetValue(dist);
            seeds->InsertElement(cnt_seeds++, node);
          }
        }
      }
    }
  }
  return {seeds, outside};
}
