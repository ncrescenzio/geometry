
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

file(GLOB vtk_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*.F90)
add_library(vtk ${vtk_SOURCES})

target_include_directories(vtk INTERFACE
   $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
   $<INSTALL_INTERFACE:mod>
)

include(GNUInstallDirs)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
   DESTINATION mod
   COMPONENT library
)
install(TARGETS vtk
   EXPORT geometry-targets
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   COMPONENT library
)
