add_definitions(-DUSE_FORTRAN_BLAS)

set (Directories
 .
 common
 common/tunicate
 common/newsparse
 LosTopos3D
)

append_files (Headers "h" ${Directories})
append_files (Sources "cpp" ${Directories})

if(NOT WIN32)
find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})
endif(NOT WIN32)

add_library(LosTopos STATIC ${Headers} ${Sources})

if(NOT WIN32)
target_link_libraries (LosTopos ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
else(NOT WIN32)
target_link_libraries (LosTopos ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif(NOT WIN32)