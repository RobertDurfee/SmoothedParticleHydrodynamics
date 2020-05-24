cmake_minimum_required(VERSION 2.8)

project(final)

find_package(OpenGL REQUIRED)

if (APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

if (UNIX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=gnu++11")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable")
  # recommended but not set by default
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
elseif(MSVC)
  # recommended but not set by default
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -WX")
endif()

set (LIBS ${OPENGL_gl_LIBRARY})

# GLFW
set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)
add_subdirectory(3rd_party/glfw)
list(APPEND LIBS glfw)
list(APPEND INCLUDES 3rd_party/glfw/include)

# GLEW - not needed on OS X
# we add glew source/header directly to the build, no glew library build.
if (NOT APPLE)
  add_definitions(-DGLEW_STATIC)
  list(APPEND INCLUDES 3rd_party/glew/include)
  list(APPEND SRC 3rd_party/glew/src/glew.c)
  SOURCE_GROUP(GLEW FILES 3rd_party/glew/src/glew.c)
endif()


# vecmath include directory
include_directories(vecmath/include)
add_subdirectory(vecmath)
list (APPEND LIBS vecmath)
list (APPEND INCLUDES vecmath/include)
list (APPEND SRC
  src/main.cpp
  src/util.cpp
  src/camera.cpp
  src/vertexrecorder.cpp
  src/glprogram.cpp
  src/sph.cpp
)
list (APPEND HEADER
  src/gl.h
  src/util.h
  src/camera.h
  src/vertexrecorder.h
  src/glprogram.h
  src/sphparticle.h
  src/sphobstacle.h
  src/sph.h
)

add_executable(final ${SRC} ${HEADER})
target_include_directories(final PUBLIC ${INCLUDES})
target_link_libraries(final ${LIBS})
