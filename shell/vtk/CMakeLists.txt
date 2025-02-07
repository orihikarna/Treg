cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

project(TrackEgg)
find_package(ITK REQUIRED)
find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersSources
  IOGeometry
  IOImage
  IOLegacy
  ImagingGeneral
  ImagingSources
  ImagingHybrid
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  RenderingVolume
  RenderingVolumeOpenGL2
)

include(${ITK_USE_FILE})
if (NOT VTK_FOUND)
  message(FATAL_ERROR "TrackEgg: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
add_executable(TrackEgg MACOSX_BUNDLE TrackEgg.cxx vtkImageEuclideanDistance2.cxx)
target_link_libraries(TrackEgg PRIVATE ${VTK_LIBRARIES})
target_link_libraries(TrackEgg PRIVATE ${ITK_LIBRARIES})
target_link_libraries(TrackEgg PRIVATE ${ITKVNL_LIBRARIES})
  
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS TrackEgg
  MODULES ${VTK_LIBRARIES}
)