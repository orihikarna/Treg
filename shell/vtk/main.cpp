#include <sys/resource.h>
#include <sys/time.h>
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

#include "TrackEgg.h"
#include "vtkImageEuclideanDistance2.h"

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

long dumpMemoryUsage(const std::string &title = "") {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    const int mem_MB = std::round(usage.ru_maxrss / 1024.0 / 1024.0);
    std::cout << "===[ Memory: " << mem_MB << " MB (" << title << ") ]===" << std::endl;
    return usage.ru_maxrss;  // bytes
  } else
    return 0;
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
    const itk::Size<3> size{SizeX, SizeY, SizeZ};
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
