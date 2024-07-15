cmake_minimum_required(VERSION 3.5.0)
project(spectrogram VERSION 0.1.0 LANGUAGES C CXX)
SET(CMAKE_BUILD_TYPE Debug)
include(CTest)
enable_testing()
SET(CMAKE_CXX_FLAGS "-lfftw3_threads -lfftw3 -lm" )
SET(CMAKE_C_FLAGS   "-lfftw3_threads -lfftw3 -lm")

INCLUDE_DIRECTORIES(/usr/include)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/source/headers)
find_package(UHD REQUIRED)
find_package(HDF5 REQUIRED)


include_directories(${UHD_INCLUDE_DIRS})

LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)
include_directories(${YOUR_DIRECTORY})
set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1) 
set(SOURCES  source/src/focus.cpp
            source/src/functions.cpp
)
add_executable(spectrogram ${SOURCES})
target_link_libraries(spectrogram ${HDF5_LIBRARIES})


target_link_libraries(spectrogram HDF5::HDF5)
target_include_directories(spectrogram PRIVATE ${HDF5_INCLUDE_DIRS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(spectrogram fftw3)


project(spectrogram VERSION 0.1.0 LANGUAGES C CXX)

SET(CMAKE_CXX_FLAGS "-lfftw3_threads -lfftw3 -lm" )
SET(CMAKE_C_FLAGS   "-lfftw3_threads -lfftw3 -lm")

INCLUDE_DIRECTORIES(/usr/include)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/source/headers)
find_package(UHD REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX) # Ensure both C and C++ libraries are found

include_directories(${UHD_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS}) # Include HDF5 headers

LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)
include_directories(${YOUR_DIRECTORY})
set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1) 
set(SOURCES  source/src/focus.cpp
            source/src/functions.cpp
)

target_link_libraries(spectrogram ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}) # Link HDF5 libraries