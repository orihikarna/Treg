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

float calc_dist2_egg(const Eigen::Vector3f &pos);
float calc_dist2_track_egg(const Eigen::Vector3f &pos);

void calc_ridge_points();

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

std::tuple<NodeContainer::Pointer, NodeContainer::Pointer> TrackEggSeeds();
