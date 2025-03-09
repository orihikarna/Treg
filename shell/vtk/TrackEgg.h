#pragma once

#include <itkFastMarchingImageFilter.h>
#include <vnl/algo/vnl_powell.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_nonlinear_minimizer.h>
#include <vtkImageData.h>

#include <Eigen/Dense>
#include <cmath>

constexpr float ball_r = 57.2 / 2;
constexpr float hole_r = ball_r + 2.5 / 2;
constexpr float btm_h = 1;

// center = ball
// center = [0, -30, hole_r + btm_h];
constexpr float ball_z = -30;
constexpr float ball_y = hole_r + btm_h;

constexpr float egg_org_y = -ball_y;
constexpr float egg_org_z = -ball_z;
const Eigen::Vector3f egg_org{0, egg_org_y, egg_org_z};

inline constexpr float rad2deg(float rad) { return rad * 180.0f / float(M_PI); }
inline constexpr float deg2rad(float deg) { return deg / 180.0f * float(M_PI); }

constexpr float egg_alpha = deg2rad(36);
constexpr float egg_tilt = deg2rad(11);
const float tilt_co = std::cos(egg_tilt);
const float tilt_si = std::sin(egg_tilt);

inline std::tuple<float, float> egg_translate_fwd(float y, float z) { return {y + egg_org_y, z + egg_org_z}; }
inline std::tuple<float, float> egg_translate_bck(float y, float z) { return {y - egg_org_y, z - egg_org_z}; }
inline std::tuple<float, float> egg_rotate_fwd(float y, float z) { return {tilt_co * y + tilt_si * z, tilt_co * z - tilt_si * y}; }
inline std::tuple<float, float> egg_rotate_bck(float y, float z) { return {tilt_co * y - tilt_si * z, tilt_co * z + tilt_si * y}; }

inline Eigen::Vector3f egg_translate_fwd(const Eigen::Vector3f &pos) {
  const auto [y, z] = egg_translate_fwd(pos[1], pos[2]);
  return Eigen::Vector3f{pos[0], y, z};
}
inline Eigen::Vector3f egg_translate_bck(const Eigen::Vector3f &pos) {
  const auto [y, z] = egg_translate_bck(pos[1], pos[2]);
  return Eigen::Vector3f{pos[0], y, z};
}
inline Eigen::Vector3f egg_rotate_fwd(const Eigen::Vector3f &pos) {
  const auto [y, z] = egg_rotate_fwd(pos[1], pos[2]);
  return Eigen::Vector3f{pos[0], y, z};
}
inline Eigen::Vector3f egg_rotate_bck(const Eigen::Vector3f &pos) {
  const auto [y, z] = egg_rotate_bck(pos[1], pos[2]);
  return Eigen::Vector3f{pos[0], y, z};
}

constexpr float mkw_r = 3.0;

constexpr float hole_mkw_r = hole_r + mkw_r;
constexpr float egg_scale_x = 38 - mkw_r;
constexpr float egg_scale_y = 42 - mkw_r;
constexpr float egg_scale_z = 66 - mkw_r;

const float egg_zbtm = 0;
const float egg_ztop = std::tan(egg_alpha);
const float egg_rad_btm = 1;
const float egg_rad_top = 2 - 1 / std::cos(egg_alpha);
const float egg_zmin = egg_zbtm - egg_rad_btm;
const float egg_zmax = egg_ztop + egg_rad_top;
const float egg_zedge_btm = 0;
const float egg_zedge_top = 2 * std::sin(egg_alpha);

constexpr float spacing = 1.0 / 1;
// constexpr float spacing = 3.0 / 4;
// constexpr float spacing = 2.0 / 3;
constexpr int SizeX = int(110 / spacing + 0.5f);
constexpr int SizeY = int(120 / spacing + 0.5f);
constexpr int SizeZ = int(220 / spacing + 0.5f);

constexpr int OrigX = SizeX / 2;
constexpr int OrigY = SizeY * 3 / 4;
constexpr int OrigZ = SizeZ / 4;

// ===== 1. inside / outside for Euclidean Distance Tranform =====
template <typename T>
void Egg(vtkImageData *img) {
  T *ptr_vol = (T *)img->GetScalarPointer();
  for (int iz = 0; iz < SizeZ; ++iz) {
    T *ptr_z = ptr_vol + SizeY * SizeX * iz;
    const float z = spacing * (iz - OrigZ);
    const float sz = z * (1 / egg_scale_z);
    const float dz = sz - egg_org_z / egg_scale_z;
    for (int iy = 0; iy < SizeY; ++iy) {
      T *ptr_y = ptr_z + SizeX * iy;
      const float y = spacing * (iy - OrigY);
      const float sy = y * (1 / egg_scale_y);
      const float sy2 = sy * sy;
      for (int ix = 0; ix < SizeX; ++ix) {
        const float x = spacing * (ix - OrigX);
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
    const float z = spacing * (iz - OrigZ);
    const float dz = z - ball_z;
    const float dz2 = dz * dz;
    for (int iy = 0; iy < SizeY; ++iy) {
      T *ptr_y = ptr_z + SizeX * iy;
      const float y = spacing * (iy - OrigY);
      const float dy = y - ball_y;
      const float dy2 = dy * dy;
      for (int ix = 0; ix < SizeX; ++ix) {
        const float x = spacing * (ix - OrigX);
        const float dx = x;
        const float dx2 = dx * dx;
        if (dx2 + dy2 + dz2 < hole_mkw_r * hole_mkw_r) {
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
