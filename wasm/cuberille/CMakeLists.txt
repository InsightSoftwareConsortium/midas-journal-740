add_executable(cuberille cuberille.cxx)
target_link_libraries(cuberille PUBLIC ${ITK_LIBRARIES})

add_test(NAME cuberille-help COMMAND cuberille --help)

add_test(NAME cuberille-defaults
  COMMAND cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/fuel.mha
    ${CMAKE_CURRENT_BINARY_DIR}/cuberille-defaults.vtk)

add_test(NAME cuberille-convergence-options
  COMMAND cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/fuel.mha
    ${CMAKE_CURRENT_BINARY_DIR}/cuberille-convergence-options.vtk
    --surface-distance-threshold 0.2
    --step-length 0.24
    --step-relaxation-factor 0.95
    --max-steps 100)

add_test(NAME cuberille-no-projection
  COMMAND cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/fuel.mha
    ${CMAKE_CURRENT_BINARY_DIR}/cuberille-no-projection.vtk
    --surface-distance-threshold 0.2
    --step-length 0.24
    --step-relaxation-factor 0.95
    --max-steps 100
    --no-projection)

add_test(NAME cuberille-quadrilateral-faces
  COMMAND cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/fuel.mha
    ${CMAKE_CURRENT_BINARY_DIR}/cuberille-quadrilateral-faces.vtk
    --surface-distance-threshold 0.2
    --step-length 0.24
    --step-relaxation-factor 0.95
    --max-steps 100
    --quadrilateral-faces)

add_test(NAME cuberille-bspline
  COMMAND cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/fuel.mha
    ${CMAKE_CURRENT_BINARY_DIR}/cuberille-bspline.vtk
    --surface-distance-threshold 0.2
    --step-length 0.24
    --step-relaxation-factor 0.95
    --max-steps 100
    --interpolator bspline)

add_test(NAME cuberille-windowed-sinc
  COMMAND cuberille
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test/Input/fuel.mha
    ${CMAKE_CURRENT_BINARY_DIR}/cuberille-windowed-sinc.vtk
    --surface-distance-threshold 0.2
    --step-length 0.24
    --step-relaxation-factor 0.95
    --max-steps 100
    --interpolator windowed-sinc)