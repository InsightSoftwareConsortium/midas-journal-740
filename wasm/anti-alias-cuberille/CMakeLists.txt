add_executable(anti-alias-cuberille anti-alias-cuberille.cxx)

set(io_components)
if (NOT WASI AND NOT EMSCRIPTEN)
  set(io_components
    ITKIONIFTI
    ITKIOMeshVTK)
endif()

find_package(ITK REQUIRED COMPONENTS
  WebAssemblyInterface
  ITKAntiAlias
  ITKSmoothing
  ITKLabelMap
  ${io_components}
)
include(${ITK_USE_FILE})

target_link_libraries(anti-alias-cuberille PUBLIC ${ITK_LIBRARIES})

add_test(NAME anti-alias-cuberille-help COMMAND anti-alias-cuberille --help)

add_test(NAME anti-alias-cuberille-defaults
  COMMAND anti-alias-cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/brain-seg.nii.gz
    ${CMAKE_CURRENT_BINARY_DIR}/anti-alias-cuberille-defaults.vtk)

add_test(NAME anti-alias-cuberille-convergence-options
  COMMAND anti-alias-cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/brain-seg.nii.gz
    ${CMAKE_CURRENT_BINARY_DIR}/anti-alias-cuberille-convergence-options.vtk
    --surface-distance-threshold 0.01
    --maximum-steps 150
    --step-relaxation-factor 0.75)

add_test(NAME anti-alias-cuberille-quadrilateral-faces
  COMMAND anti-alias-cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/brain-seg.nii.gz
    ${CMAKE_CURRENT_BINARY_DIR}/anti-alias-cuberille-quadrilateral-faces.vtk
    --quadrilateral-faces
    --surface-distance-threshold 0.5
    --maximum-steps 150
    --step-relaxation-factor 0.75
    )

add_test(NAME anti-alias-cuberille-bspline
  COMMAND anti-alias-cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/brain-seg.nii.gz
    ${CMAKE_CURRENT_BINARY_DIR}/anti-alias-cuberille-bspline.vtk
    --surface-distance-threshold 0.5
    --maximum-steps 200
    --interpolator bspline)

add_test(NAME anti-alias-cuberille-windowed-sinc
  COMMAND anti-alias-cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/brain-seg.nii.gz
    ${CMAKE_CURRENT_BINARY_DIR}/anti-alias-cuberille-windowed-sinc.vtk
    --surface-distance-threshold 0.5
    --maximum-steps 200
    --interpolator windowed-sinc)