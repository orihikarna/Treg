#include <itkFastMarchingImageFilter.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vnl/algo/vnl_powell.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_nonlinear_minimizer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkImageData.h>
#include <vtkImageEuclideanDistance.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageShiftScale.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSTLWriter.h>
#include <vtkVersion.h>
#include <vtkVolumeProperty.h>
#include <vtkVoxelModeller.h>

#include <Eigen/Dense>
#include <sstream>

#include "vtkImageEuclideanDistance2.h"

long dumpMemoryUsage(const std::string &title = "") {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    const int mem_MB = std::round(usage.ru_maxrss / 1024.0 / 1024.0);
    std::cout << "===[ Memory: " << mem_MB << " MB (" << title << ") ]===" << std::endl;
    return usage.ru_maxrss;  // bytes
  } else
    return 0;
}

// vtkFlyingEdges3D was introduced in VTK >= 8.2
#if VTK_MAJOR_VERSION >= 9 || (VTK_MAJOR_VERSION >= 8 && VTK_MINOR_VERSION >= 2)
#define USE_FLYING_EDGES
#else
#undef USE_FLYING_EDGES
#endif

#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
#else
#include <vtkMarchingCubes.h>
#endif

constexpr unsigned int Dimension = 3;
using InternalPixelType = float;
using InternalImageType = itk::Image<InternalPixelType, Dimension>;
using FastMarchingFilterType = itk::FastMarchingImageFilter<InternalImageType, InternalImageType>;
using NodeContainer = FastMarchingFilterType::NodeContainer;
using NodeType = FastMarchingFilterType::NodeType;

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
        T val = -1;
        // val = -dist;
        val = (dist < 0) ? 0 : -1;
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
                // dist = std::max(fine_dist_egg_pxl, dist_ball_pxl);
                // dist = std::sqrt(fine_dist_egg_pxl * fine_dist_egg_pxl + dist_ball_pxl * dist_ball_pxl);
              }
            }
            // if (-fine_dist_egg_pxl < dist_ball_pxl) {  // use egg dist
            //   node.SetValue(fine_dist_egg_pxl;
            // } else {  // use ball dist
            //   node.SetValue(-dist_ball_pxl);
            // }
            node.SetValue(dist);
            seeds->InsertElement(cnt_seeds++, node);
          }
        }
      }
    }
  }
  return {seeds, outside};
}

vtkNew<vtkImageData> createImageData(int scalar_type) {
  vtkNew<vtkImageData> data;
  data->SetDimensions(SizeX, SizeY, SizeZ);
  data->SetSpacing(spacing, spacing, spacing);
  data->SetOrigin(-SizeX / 2 * spacing, -SizeY / 2 * spacing, -SizeZ / 2 * spacing);
  data->AllocateScalars(scalar_type, 1);
  const size_t mem_size = size_t(data->GetScalarSize()) * SizeX * SizeY * SizeZ;
  const size_t mem_size_MB = mem_size / 1024 / 1024;
  std::cout << SizeX << ", " << SizeY << ", " << SizeZ << std::endl;
  std::cout << "created image data: " << mem_size_MB << " MB" << std::endl;
  dumpMemoryUsage("created image data");
  return data;
}

int main(int argc, char *argv[]) {
  vtkNew<vtkNamedColors> colors;

  dumpMemoryUsage("main");
  vtkNew<vtkImageData> dist_data_filter;
  dumpMemoryUsage("main");
  if (false) {
    vtkNew<vtkImageData> dist_data_sht;
    {
      vtkNew<vtkImageData> dist_data_flt = createImageData(VTK_FLOAT);
      {
        //   std::cout << "egg voxels" << std::endl;
        vtkNew<vtkImageData> imgdata = createImageData(VTK_CHAR);
        Egg<char>(imgdata);
        Ball<char>(imgdata);
        dumpMemoryUsage("set egg");

        std::cout << "distance" << std::endl;
        vtkNew<vtkImageEuclideanDistance2> distance;
        distance->SetInitialize(true);
        distance->SetAlgorithmToSaito();
        distance->SetInputData(imgdata);
        distance->SetConsiderAnisotropy(false);
        distance->SetOutput(dist_data_flt.Get());
        distance->Update();
        dumpMemoryUsage("distance");
      }
      dumpMemoryUsage("distance out");

      dist_data_sht = createImageData(VTK_SHORT);
      dumpMemoryUsage("distance short");
      {
        const size_t size = size_t(SizeX * SizeY) * size_t(SizeZ);
        const DistType *pin = (DistType *)dist_data_flt->GetScalarPointer();
        short *pout = (short *)dist_data_sht->GetScalarPointer();
        for (size_t n = 0; n < size; ++n) {
          pout[n] = -short(std::round(std::min(std::sqrt(pin[n]) * 1024, DistType(32767))));
        }
      }
      dumpMemoryUsage("sqrt");

      if (false) {
        typedef short T;
        T *ptr_vol = (T *)dist_data_sht->GetScalarPointer();
        const int iz = SizeZ / 2;
        T *ptr_z = ptr_vol + SizeY * SizeX * iz;
        const int iy = SizeY / 2;
        T *ptr_y = ptr_z + SizeX * iy;
        for (int ix = 0; ix < SizeX; ++ix) {
          const T val = ptr_y[ix];
          if (val) {
            std::cout << "ix = " << ix << ", dist = " << std::sqrt(val) << std::endl;
          }
        }
      }
      if (false) {
        typedef short T;
        T *ptr_vol = (T *)dist_data_sht->GetScalarPointer();
        const int iz = SizeZ / 2;
        T *ptr_z = ptr_vol + SizeY * SizeX * iz;
        for (int iy = 0; iy < SizeY; ++iy) {
          T *ptr_y = ptr_z + SizeX * iy;
          const int ix = SizeX / 2;
          const T val = ptr_y[ix];
          if (val) {
            std::cout << "iy = " << iy << ", dist = " << std::sqrt(val) << std::endl;
          }
        }
      }
    }

    dumpMemoryUsage("distance short out");

    std::cout << "gaussian smoothing" << std::endl;

    dist_data_filter = createImageData(VTK_SHORT);
    {
      vtkNew<vtkImageGaussianSmooth> smoothing_filter;
      smoothing_filter->SetDimensionality(3);
      smoothing_filter->SetInputData(dist_data_sht);
      constexpr double sigma = 3.0;
      constexpr double radius = 4;
      smoothing_filter->SetStandardDeviations(sigma, sigma, sigma);
      smoothing_filter->SetRadiusFactors(radius, radius, radius);
      smoothing_filter->SetOutput(dist_data_filter);
      smoothing_filter->Update();
      dumpMemoryUsage("smoothed");
    }
  } else {
    dumpMemoryUsage("trial points");
    auto fastMarching = FastMarchingFilterType::New();
    if (false) {
      auto seeds = NodeContainer::New();
      seeds->Initialize();
      constexpr double M = 1.5;
      int cnt = 0;
      for (double dz = -M; dz <= M; dz += 1) {
        for (double dy = -M; dy <= M; dy += 1) {
          for (double dx = -M; dx <= M; dx += 1) {
            InternalImageType::IndexType seedPosition;
            seedPosition[0] = SizeX / 2 + dx;
            seedPosition[1] = SizeY / 2 + dy;
            seedPosition[2] = SizeZ / 2 + dz;
            NodeType node;
            const double seedValue = sqrt(dx * dx + dy * dy + dz * dz);
            node.SetValue(seedValue);
            node.SetIndex(seedPosition);
            seeds->InsertElement(cnt++, node);
          }
        }
      }
      fastMarching->SetTrialPoints(seeds);
    } else {
      auto [seeds, outside] = TrackEggSeeds();
      fastMarching->SetTrialPoints(seeds);
      fastMarching->SetOutsidePoints(outside);
    }
    dumpMemoryUsage("fast marching");
    const double stoppingTime = 2 * mkw_r / spacing;
    const itk::Size<Dimension> size{SizeX, SizeY, SizeZ};
    fastMarching->SetOutputSize(size);
    fastMarching->SetSpeedConstant(1.0);
    fastMarching->SetStoppingValue(stoppingTime);
    fastMarching->Update();

    dumpMemoryUsage("convert distance");
    auto dist = fastMarching->GetOutput();
    dist_data_filter = createImageData(VTK_FLOAT);
    const size_t num = size_t(SizeX * SizeY) * SizeZ;
    const float *pin = dist->GetBufferPointer();
    float *pout = (float *)dist_data_filter->GetScalarPointer();
    for (size_t n = 0; n < num; ++n) {
      pout[n] = -pin[n] / spacing * 1024;
      // pout[n] = -pin[n];
    }
  }

  dumpMemoryUsage("smoothed out");

  vtkNew<vtkRenderer> renderer;
  if (false) {
    std::cout << "volume rendering" << std::endl;

    vtkNew<vtkPiecewiseFunction> opacityTransferFunction;
    // Create transfer mapping scalar value to opacity
    if (false) {
      opacityTransferFunction->AddPoint(-1, 0);
      opacityTransferFunction->AddPoint(0, 1);
    } else {
      opacityTransferFunction->AddPoint(2, 1);
      opacityTransferFunction->AddPoint(3, 0);
    }

    // Create transfer mapping scalar value to color
    vtkNew<vtkColorTransferFunction> colorTransferFunction;
    if (true) {
      colorTransferFunction->AddRGBPoint(-32768, 0.0, 1, 1);
      colorTransferFunction->AddRGBPoint(+32767, 0.0, 1, 1);
    } else {
      colorTransferFunction->AddRGBPoint(-1024, 1.0, 0.0, 0.0);
      float div = 128;
      colorTransferFunction->AddRGBPoint(div * 0, 1.0, 0.0, 0.0);
      colorTransferFunction->AddRGBPoint(div * 1, 0.5, 0.5, 0.0);
      colorTransferFunction->AddRGBPoint(div * 2, 0.0, 1.0, 0.0);
      colorTransferFunction->AddRGBPoint(div * 3, 0.0, 0.5, 0.5);
      colorTransferFunction->AddRGBPoint(div * 4, 0.0, 0.0, 1.0);
      colorTransferFunction->AddRGBPoint(div * 5, 0.5, 0.0, 0.5);
      colorTransferFunction->AddRGBPoint(div * 6, 0.5, 0.5, 0.5);
    }

    // The property describes how the data will look
    vtkNew<vtkVolumeProperty> volumeProperty;
    volumeProperty->SetColor(colorTransferFunction);
    volumeProperty->SetScalarOpacity(opacityTransferFunction);
    volumeProperty->ShadeOn();
    volumeProperty->SetInterpolationTypeToLinear();

    // The mapper / ray cast function know how to render the data
    vtkNew<vtkFixedPointVolumeRayCastMapper> volumeMapper;
    volumeMapper->SetInputData(dist_data_filter);
    // volumeMapper->SetInputData(imgdata);

    // The volume holds the mapper and the property and
    // can be used to position/orient the volume
    vtkNew<vtkVolume> volume;
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);

    renderer->AddVolume(volume);
    renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
    // renderer->GetActiveCamera()->Azimuth(45);
    // renderer->GetActiveCamera()->Elevation(30);
    // renderer->ResetCameraClippingRange();
    // renderer->ResetCamera();
  } else {
    std::cout << "marching cube" << std::endl;
#ifdef USE_FLYING_EDGES
    vtkNew<vtkFlyingEdges3D> surface;
#else
    vtkNew<vtkMarchingCubes> surface;
#endif
    const double dist_voxel = mkw_r / spacing * 1024;
    const double isoValue = -dist_voxel;
    // const double isoValue = 0;
    surface->SetInputData(dist_data_filter);
    surface->ComputeNormalsOn();
    surface->SetValue(0, isoValue);
    dumpMemoryUsage("marching cubed");

    std::cout << "STL export" << std::endl;
    vtkNew<vtkSTLWriter> stl_writer;
    stl_writer->SetInputConnection(surface->GetOutputPort());
    std::stringstream ss;
    ss << "surface-mkw=" << mkw_r << ".stl";
    stl_writer->SetFileName(ss.str().c_str());
    stl_writer->Write();
    dumpMemoryUsage("wrote STL");

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(surface->GetOutputPort());
    mapper->ScalarVisibilityOff();

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

    renderer->AddActor(actor);
    renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());
  }
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(2000, 2000);
  renderWindow->SetWindowName("SimpleRayCast");
  renderWindow->Render();

  vtkNew<vtkRenderWindowInteractor> interactor;
  interactor->SetRenderWindow(renderWindow);
  interactor->Start();

  return EXIT_SUCCESS;
}
