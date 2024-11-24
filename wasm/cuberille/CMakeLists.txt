add_executable(cuberille cuberille.cxx)
target_link_libraries(cuberille PUBLIC ${ITK_LIBRARIES})

add_test(NAME cuberille-help COMMAND cuberille --help)
