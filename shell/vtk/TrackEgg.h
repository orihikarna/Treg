#pragma once

#include <itkFastMarchingImageFilter.h>
#include <vnl/algo/vnl_powell.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_nonlinear_minimizer.h>
#include <vtkImageData.h>

#include <Eigen/Dense>
#include <cmath>

constexpr float mkw_r = 3.0;
constexpr float btm_h = 1;

constexpr float egg_org_z = -25;

constexpr float ball_r = 57.2 / 2;
constexpr float hole_r = ball_r + 2.5 / 2 + mkw_r;

constexpr float ball_z = 80 + egg_org_z;
constexpr float ball_y = ball_r + btm_h;

constexpr float egg_scale_x = 47 - mkw_r;
constexpr float egg_scale_y = 50 - mkw_r;
constexpr float egg_scale_z = 80 - mkw_r;

// constexpr float egg_scale_x = 47 - mkw_r;
// constexpr float egg_scale_y = 50 - mkw_r;
// constexpr float egg_scale_z = 80 - mkw_r;

// constexpr float egg_alpha = 36 / 180.0f * M_PI;
constexpr float egg_alpha = 42 / 180.0f * M_PI;

const float egg_zbtm = 0;
const float egg_ztop = std::tan(egg_alpha);
const float egg_rad_btm = 1;
const float egg_rad_top = 2 - 1 / std::cos(egg_alpha);
const float egg_zmin = egg_zbtm - egg_rad_btm;
const float egg_zmax = egg_ztop + egg_rad_top;
const float egg_zedge_btm = 0;
const float egg_zedge_top = 2 * std::sin(egg_alpha);

constexpr float spacing = 1.0 / 1;
// constexpr float spacing = 2.0 / 3;
constexpr int SizeX = int(110 / spacing + 0.5f);
constexpr int SizeY = int(120 / spacing + 0.5f);
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

// ridge points
// center = ball center
// radius = ball radius
// azim = -PI ~ +PI
class TrackEggRidge {
 private:
  const size_t M;
  const size_t N;

  using T = float;
  const T dazim = T(M_PI) / M;

  std::vector<T> azims_;
  std::vector<T> elevs_;
  std::vector<Eigen::Vector3f> points_;

 public:
  TrackEggRidge(size_t _M = 18) : M(_M), N(2 * M + 1), dazim(T(M_PI) / M) { Init(); }

  bool IsOnHoleSurface(float azim, float elev) const;
  T CalcMinDist(const Eigen::Vector3f &pos) const;

 private:
  void Init();
};

using InternalImageType = itk::Image<float, 3>;
using FastMarchingFilterType = itk::FastMarchingImageFilter<InternalImageType, InternalImageType>;
using NodeContainer = FastMarchingFilterType::NodeContainer;
using NodeType = FastMarchingFilterType::NodeType;

std::tuple<NodeContainer::Pointer, NodeContainer::Pointer> TrackEggSeeds();
