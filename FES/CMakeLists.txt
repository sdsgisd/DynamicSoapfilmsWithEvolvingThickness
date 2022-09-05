append_files (Headers "h" .)
append_files (Sources "cpp" .)

add_executable (FES ${Headers} ${Sources})

TARGET_LINK_LIBRARIES(FES LosTopos ${GLEW_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
