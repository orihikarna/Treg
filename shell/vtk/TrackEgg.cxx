#include <itkFastMarchingImageFilter.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkDICOMImageReader.h>
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
#include <vtkSphereSource.h>
#include <vtkVersion.h>
#include <vtkVolumeProperty.h>
#include <vtkVoxelModeller.h>

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

constexpr float mkw_r = 3.0;
constexpr float btm_h = 1;

constexpr float egg_origin_z = -25;

constexpr float ball_r = 57.2 / 2;
constexpr float hole_r = ball_r + 2.5 / 2 + mkw_r;

constexpr float ball_z = 80 + egg_origin_z;
constexpr float ball_y = ball_r + btm_h;

constexpr float egg_scale_x = 45 - mkw_r;
constexpr float egg_scale_y = 40 - mkw_r;
constexpr float egg_scale_z = 80 - mkw_r;

constexpr float egg_alpha = 36 / 180.0f * M_PI;
const float egg_org_top = std::tan(egg_alpha);
const float egg_rad_top = 2 - 1 / std::cos(egg_alpha);
const float egg_zh = 2 * std::sin(egg_alpha);

constexpr float spacing = 1.0 / 1;
constexpr int SizeX = int(110 / spacing + 0.5f);
constexpr int SizeY = int(100 / spacing + 0.5f);
constexpr int SizeZ = int(220 / spacing + 0.5f);

float egg_slice_r2(float z) {
  if (z < -0.5f) return -1;
  if (z < 0) return 0.25f - z * z;
  if (z < egg_zh) {
    const float d = std::sqrt(1 - z * z) - 0.5f;
    return d * d;
  }
  const float dz = z - egg_org_top;
  if (dz < egg_rad_top) return egg_rad_top * egg_rad_top - dz * dz;
  return -1;
}

void Egg_slice(vtkImageData *img) {
  int16_t *ptr_vol = (int16_t *)img->GetScalarPointer();
  for (int iz = 0; iz < SizeZ; ++iz) {
    int16_t *ptr_z = ptr_vol + SizeY * SizeX * iz;
    const float z = spacing * (iz - SizeZ / 2);
    const float r2 = egg_slice_r2((z - egg_origin_z / spacing) / egg_scale_z);
    for (int iy = 0; iy < SizeY; ++iy) {
      int16_t *ptr_y = ptr_z + SizeX * iy;
      const float y = spacing * (iy - SizeY / 2);
      const float sy = y * (1 / egg_scale_y);
      const float sy2 = sy * sy;
      for (int ix = 0; ix < SizeX; ++ix) {
        const float x = spacing * (ix - SizeX / 2);
        const float sx = x * (1 / egg_scale_x);
        const float sx2 = sx * sx;
        ptr_y[ix] = std::min<float>(std::max<float>((sx2 + sy2 - r2) * 1024, -16), +16);
      }
    }
  }
}

template <typename T>
void Egg(vtkImageData *img) {
  T *ptr_vol = (T *)img->GetScalarPointer();
  for (int iz = 0; iz < SizeZ; ++iz) {
    T *ptr_z = ptr_vol + SizeY * SizeX * iz;
    const float z = spacing * (iz - SizeZ / 2);
    const float sz = z * (1 / egg_scale_z);
    const float dz = sz - egg_origin_z / egg_scale_z;
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
            const float dz_top = dz - egg_org_top;
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
  if (true) {
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
    constexpr unsigned int Dimension = 3;
    using InternalPixelType = float;
    using InternalImageType = itk::Image<InternalPixelType, Dimension>;
    using FastMarchingFilterType = itk::FastMarchingImageFilter<InternalImageType, InternalImageType>;
    auto fastMarching = FastMarchingFilterType::New();

    using NodeContainer = FastMarchingFilterType::NodeContainer;
    using NodeType = FastMarchingFilterType::NodeType;
    {
      auto seeds = NodeContainer::New();
      InternalImageType::IndexType seedPosition;
      seedPosition[0] = SizeX / 2;
      seedPosition[1] = SizeY / 2;
      seedPosition[2] = SizeZ / 2;
      NodeType node;
      constexpr double seedValue = 0.0;
      node.SetValue(seedValue);
      node.SetIndex(seedPosition);
      seeds->Initialize();
      seeds->InsertElement(0, node);
      fastMarching->SetTrialPoints(seeds);
    }
    const itk::Size<Dimension> size{SizeX, SizeY, SizeZ};
    fastMarching->SetOutputSize(size);
    const double stoppingTime = mkw_r * 2;
    fastMarching->SetSpeedConstant(1.0);
    fastMarching->SetStoppingValue(stoppingTime);
    fastMarching->Update();
    // auto mapWriter = InternalWriterType::New();
    // mapWriter->SetInput(fastMarching->GetOutput());
    // mapWriter->SetFileName("FastMarchingFilterOutput4.mha");
    // mapWriter->Update();

    auto dist = fastMarching->GetOutput();
    dist_data_filter = createImageData(VTK_FLOAT);
    const size_t num = size_t(SizeX * SizeY) * SizeZ;
    const float *pin = dist->GetBufferPointer();
    float *pout = (float *)dist_data_filter->GetScalarPointer();
    for (size_t n = 0; n < num; ++n) {
      pout[n] = -pin[n] / spacing * 1024;
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
