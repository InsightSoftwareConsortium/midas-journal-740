cmake_minimum_required(VERSION 3.16)
project(itkwasm-cuberille LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(io_components)
if (NOT WASI AND NOT EMSCRIPTEN)
  set(io_components
    ITKIOMeta
    ITKIOMeshVTK)
endif()

find_package(ITK REQUIRED COMPONENTS
  WebAssemblyInterface
  ${io_components}
)
include(${ITK_USE_FILE})

enable_testing()

# Begin create-itk-wasm added pipelines.
add_subdirectory(cuberille)
add_subdirectory(anti-alias-cuberille)
# End create-itk-wasm added pipelines.
