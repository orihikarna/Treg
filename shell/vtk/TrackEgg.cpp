#include "TrackEgg.h"

namespace {
const Eigen::Vector3f ball_ctr{0, ball_y, ball_z};
const Eigen::Vector3f EGG_OOB = Eigen::Vector3f::Zero();
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

class DistanceToTrackEggSurface : public vnl_cost_function {
 private:
  const Eigen::Vector3f ball_ctr_{0, ball_y, ball_z};
  Eigen::Vector3f pos_;

 public:
  DistanceToTrackEggSurface(const Eigen::Vector3f &pos) : pos_(pos), vnl_cost_function(2) {}

  double f(vnl_vector<double> const &x) override {
    const float z = x[0];
    const float theta = x[1];
    const Eigen::Vector3f pos_egg = egg_surface(z, theta);
    if (pos_egg == EGG_OOB) return 1e9;
    const float dist2 = (pos_ - pos_egg).squaredNorm();
    const float dist_ball = hole_r - (pos_egg - ball_ctr_).norm();
    return dist2 + 1e2 * dist_ball * dist_ball;
  }
};

// ===== 2. distance to egg surface for Fast Marching =====
float calc_dist2_egg(const Eigen::Vector3f &pos) {
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
  const float dist2 = minimizer.get_end_error();
  return dist2;
}

float calc_dist2_track_egg(const Eigen::Vector3f &pos) {
  DistanceToTrackEggSurface dist(pos);
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
  // const float dist2 = minimizer.get_end_error();
  const float dist2 = (egg_surface(x[0], x[1]) - pos).squaredNorm();
  return dist2;
}

Eigen::Vector3f calc_ball_pos(float azim, float elev, float r = hole_r) {
  Eigen::Vector3f pos;
  pos[2] = (r * std::cos(elev)) * std::cos(azim);
  pos[0] = (r * std::cos(elev)) * std::sin(azim);
  pos[1] = (r * std::sin(elev));
  return pos + ball_ctr;
}

class DistanceToEggSurface_BallAzim_1pass : public vnl_cost_function {
 private:
  float azim_;

 public:
  DistanceToEggSurface_BallAzim_1pass(float azim) : azim_(azim), vnl_cost_function(3) {}

  double f(vnl_vector<double> const &x) override {
    const float z = x[0];
    const float theta = x[1];
    const float elev = x[2];
    if (elev < -M_PI / 2 || +M_PI / 2 < elev) return 1e9;
    const Eigen::Vector3f pos_egg = egg_surface(z, theta);
    if (pos_egg == EGG_OOB) return 1e9;
    const Eigen::Vector3f pos_ball = calc_ball_pos(azim_, elev);
    const float dist2 = (pos_ball - pos_egg).squaredNorm();
    return dist2;
  }
};

float calc_ball_elev_1pass(float azim) {
  DistanceToEggSurface_BallAzim_1pass dist(azim);
  vnl_powell minimizer(&dist);
  minimizer.set_f_tolerance(1e-8);
  minimizer.set_trace(true);
  const float elev0 = 60 * (M_PI / 180);
  const Eigen::Vector3f pos_ball0 = calc_ball_pos(azim, elev0);
  vnl_vector<double> x(3);
  x[0] = pos_ball0[2];  // z
  x[1] = 0;             // theta
  x[2] = elev0;
  minimizer.minimize(x);
  const float dist2 = minimizer.get_end_error();
  if (dist2 > 1e-5) {
    std::cerr << "[ERROR] dist2 = " << dist2 << std::endl;
  }
  return x[2];
}

class DistanceToEggSurface_BallAzim_2pass : public vnl_cost_function {
 private:
  float azim_;

 public:
  DistanceToEggSurface_BallAzim_2pass(float azim) : azim_(azim), vnl_cost_function(1) {}

  double f(vnl_vector<double> const &x) override {
    const float elev = x[0];
    if (elev < -M_PI / 2 || +M_PI / 2 < elev) return 1e9;
    const Eigen::Vector3f pos_ball = calc_ball_pos(azim_, elev);
    Eigen::Vector3f pos_egg;
    {
      DistanceToEggSurface dist(pos_ball);
      vnl_powell minimizer(&dist);
      minimizer.set_f_tolerance(1e-6);
      minimizer.set_trace(true);
      double z = pos_ball[2];
      if (false) {  // bring it inside the egg
        z = std::max(z, egg_zmin + egg_org_z + 1.0);
        z = std::min(z, egg_zmax + egg_org_z - 1.0);
      }
      vnl_vector<double> _x(2);
      _x[0] = z;
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
  const float elev0 = -30 * (M_PI / 180);
  vnl_vector<double> x(1);
  x[0] = elev0;
  minimizer.minimize(x);
  const float dist2 = minimizer.get_end_error();
  if (dist2 > 1e-5) {
    std::cerr << "[ERROR] dist2 = " << dist2 << std::endl;
  }
  return x[0];
}

void TrackEggRidge::Init() {
  azims_.resize(N);
  elevs_.resize(N);
  for (size_t n = 0; n < N; ++n) {
    const int m = n - M;
    const T azim = dazim * m;
    const T elev = calc_ball_elev_2pass(azim);
    azims_[n] = azim;
    elevs_[n] = elev;
    // std::cout << "(azim, elev) = " << azim * (180 / M_PI) << ", " << elev * (180 / M_PI) << std::endl;
  }
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
      const float by2 = by * by;               // baaa^2
      for (int ix = 0; ix < SizeX; ++ix) {
        pos[0] = ix;
        const float x = spacing * (ix - SizeX / 2);
        const float ex = x * (1 / egg_scale_x);  // scaled for egg
        const float ex2 = ex * ex;               // egg^2
        const float bx = x;                      // ball
        const float bx2 = bx * bx;               // ball^2

        // egg
        const float r2 = ex2 + ey2;
        float dist_egg = std::numeric_limits<float>::max();
        if (ez <= 0) {
          dist_egg = std::sqrt(r2 + ez2) - 1;
        } else {
          const float R = std::sqrt(r2) + 1;
          if (std::atan2(ez, R) < egg_alpha) {
            dist_egg = std::sqrt(R * R + ez2) - 2;
          } else {
            const float dz_top = ez - egg_ztop;
            dist_egg = std::sqrt(r2 + dz_top * dz_top) - egg_rad_top;
          }
        }
        const float dist_egg_pxl = dist_egg * min_scale / spacing;  // possible min distance
        // inside egg negative, outside egg positive

        // ball
        const float dist_ball = hole_r - std::sqrt(bx2 + by2 + bz2);  // inside positive, outside negative
        const float dist_ball_pxl = dist_ball / spacing;

        if (dist_egg_pxl < 2 && dist_ball_pxl < 2) {
          NodeType node;
          node.SetIndex(pos);
          if (dist_egg_pxl < -3.6 && dist_ball_pxl < -3.6) {
            // do nothing
          } else if (dist_egg_pxl < -1.8 && dist_ball_pxl < -1.8) {
            // outside
            outside->InsertElement(cnt_outside++, node);
          } else {
            float fine_dist_egg_pxl = (dist_egg_pxl >= 0) ? +1 : -1;
            fine_dist_egg_pxl *= std::sqrt(calc_dist2_egg({x, y, z})) / spacing;
            const float dist_edge_pxl = std::sqrt(calc_dist2_track_egg({x, y, z})) / spacing;
            float dist;
            if (dist_ball_pxl < 0) {        // outside the ball
              if (fine_dist_egg_pxl < 0) {  // inside the egg
                dist = std::max(fine_dist_egg_pxl, dist_ball_pxl);
              } else {  // outside the egg
                dist = fine_dist_egg_pxl;
              }
            } else {                        // inside the ball
              if (fine_dist_egg_pxl < 0) {  // inside the egg
                dist = dist_ball_pxl;
              } else {  // outside the egg
                dist = dist_edge_pxl;
              }
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
