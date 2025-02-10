#pragma once

#include <itkFastMarchingImageFilter.h>
#include <vnl/algo/vnl_powell.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_nonlinear_minimizer.h>
#include <vtkImageData.h>

#include <Eigen/Dense>
#include <cmath>

constexpr float mkw_r = 2.0;
constexpr float btm_h = 1;

constexpr float egg_org_z = -25;

constexpr float ball_r = 57.2 / 2;
constexpr float hole_r = ball_r + 2.5 / 2 + mkw_r;

constexpr float ball_z = 80 + egg_org_z;
constexpr float ball_y = ball_r + btm_h;

constexpr float egg_scale_x = 45 - mkw_r;
constexpr float egg_scale_y = 40 - mkw_r;
constexpr float egg_scale_z = 80 - mkw_r;

constexpr float egg_alpha = 36 / 180.0f * M_PI;

const float egg_zbtm = 0;
const float egg_ztop = std::tan(egg_alpha);
const float egg_rad_btm = 1;
const float egg_rad_top = 2 - 1 / std::cos(egg_alpha);
const float egg_zmin = egg_zbtm - egg_rad_btm;
const float egg_zmax = egg_ztop + egg_rad_top;
const float egg_zedge_btm = 0;
const float egg_zedge_top = 2 * std::sin(egg_alpha);

// constexpr float spacing = 1.0 / 2;
constexpr float spacing = 2.0 / 3;
constexpr int SizeX = int(110 / spacing + 0.5f);
constexpr int SizeY = int(100 / spacing + 0.5f);
constexpr int SizeZ = int(220 / spacing + 0.5f);

// ===== 1. inside / outside for Euclidean Distance Tranform =====
template <typename T>
void Egg(vtkImageData *img) {
  T *ptr_vol = (T *)img->GetScalarPointer();
  for (int iz = 0; iz < SizeZ; ++iz) {
    T *ptr_z = ptr_vol + SizeY * SizeX * iz;
    const float z = spacing * (iz - SizeZ / 2);
    const float sz = z * (1 / egg_scale_z);
    const float dz = sz - egg_org_z / egg_scale_z;
    for (int iy = 0; iy < SizeY; ++iy) {
      T *ptr_y = ptr_z + SizeX * iy;
      const float y = spacing * (iy - SizeY / 2);
      const float sy = y * (1 / egg_scale_y);
      const float sy2 = sy * sy;
      for (int ix = 0; ix < SizeX; ++ix) {
        const float x = spacing * (ix - SizeX / 2);
        const float sx = x * (1 / egg_scale_x);
        const float sx2 = sx * sx;
        const float r2 = sx2 + sy2;
        float dist = std::numeric_limits<float>::max();
        if (dz <= 0) {
          dist = std::sqrt(r2 + dz * dz) - 1;
        } else {
          const float R = std::sqrt(r2) + 1;
          if (std::atan2(dz, R) < egg_alpha) {
            dist = std::sqrt(R * R + dz * dz) - 2;
          } else {
            const float dz_top = dz - egg_ztop;
            dist = std::sqrt(r2 + dz_top * dz_top) - egg_rad_top;
          }
        }
        const T val = (dist < 0) ? 0 : -1;
        ptr_y[ix] = val;
      }
    }
  }
}

template <typename T>
void Ball(vtkImageData *img) {
  T *ptr_vol = (T *)img->GetScalarPointer();
  for (int iz = 0; iz < SizeZ; ++iz) {
    T *ptr_z = ptr_vol + SizeY * SizeX * iz;
    const float z = spacing * (iz - SizeZ / 2);
    const float dz = z - ball_z;
    const float dz2 = dz * dz;
    for (int iy = 0; iy < SizeY; ++iy) {
      T *ptr_y = ptr_z + SizeX * iy;
      const float y = spacing * (iy - SizeY / 2);
      const float dy = y - ball_y;
      const float dy2 = dy * dy;
      for (int ix = 0; ix < SizeX; ++ix) {
        const float x = spacing * (ix - SizeX / 2);
        const float dx = x;
        const float dx2 = dx * dx;
        if (dx2 + dy2 + dz2 < hole_r * hole_r) {
          ptr_y[ix] = -1;
        }
      }
    }
  }
}

// ===== 2. distance to egg surface for Fast Marching =====
const Eigen::Vector3f EGG_OOB = Eigen::Vector3f::Zero();

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

/*

0. 稜線を求める
1. egg_shape への最短位置を求める ball の外なら有効（最短距離＝ egg 面に垂直）
2. ball への最短位置を求める。稜線より下なら有効（最短距離＝ ball 面に垂直）
3-a. 1か2が両方も有効なら短い方を採用
3-b. 1か2の片方が有効ならそれを採用？
3-c, いずれも無効なら、稜線への最短位置を求める

片方が有効の時に、稜線を確認する必要があるか？
egg 面の外側は凸麺なのでない
ball 側の外側も凸麺なのでない

*/

using InternalImageType = itk::Image<float, 3>;
using FastMarchingFilterType = itk::FastMarchingImageFilter<InternalImageType, InternalImageType>;
using NodeContainer = FastMarchingFilterType::NodeContainer;
using NodeType = FastMarchingFilterType::NodeType;

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
