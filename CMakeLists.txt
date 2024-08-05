cmake_minimum_required(VERSION 3.5.0)
project(spectrogram VERSION 0.1.0 LANGUAGES C CXX)
SET(CMAKE_BUILD_TYPE Debug)

# #############################################
# Gtest
# #############################################
set(BUILD_TESTS off)
include(GoogleTest)
include(CTest)

enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

if(BUILD_TESTS)
    message("Tests are enabled")
    include(GoogleTest)
include(CTest)

enable_testing()
    # Add the test executable
    SET(CMAKE_CXX_FLAGS "-lfftw3_threads -lfftw3 -lm")
    SET(CMAKE_C_FLAGS "-lfftw3_threads -lfftw3 -lm")

    INCLUDE_DIRECTORIES(/usr/include)
    include_directories(${HDF5_INCLUDE_DIRS})

    include_directories(${CMAKE_SOURCE_DIR}/source/headers)

    find_package(Boost REQUIRED COMPONENTS  iostreams system filesystem thread program_options)
    find_package(UHD REQUIRED)
    find_package(HDF5 REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS C CXX) # Ensure both C and C++ libraries are found
    include_directories(${BOOST_INCLUDE_DIRS})
    include_directories(${UHD_INCLUDE_DIRS})

    LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)
    include_directories(${YOUR_DIRECTORY})
    set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1)
    set(SOURCES source/test/test.cpp
                source/src/functions.cpp
    )
    
    add_executable(spectrogram ${SOURCES})
    target_link_libraries(spectrogram ${HDF5_LIBRARIES})

    target_link_libraries(spectrogram HDF5::HDF5)
    target_include_directories(spectrogram PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(spectrogram ${Boost_LIBRARIES})
    set(CPACK_PROJECT_NAME ${PROJECT_NAME})
    set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
    include(CPack)
    find_package(GTest REQUIRED)
    target_link_libraries(spectrogram fftw3)
    include_directories(${GTEST_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(/usr/include)
    include_directories(${HDF5_INCLUDE_DIRS})
    project(spectrogram VERSION 0.1.0 LANGUAGES C CXX)
    target_link_libraries(spectrogram GTest::GTest GTest::Main)
    add_test(NAME ifftTesting COMMAND spectrogram)
    add_test(NAME 5gTest COMMAND spectrogram)
    target_link_libraries(spectrogram ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}) # Link HDF5 libraries
else()
    message("Tests are disabled")

    # ##############################################
    SET(CMAKE_CXX_FLAGS "-lfftw3_threads -lfftw3 -lm")
    SET(CMAKE_C_FLAGS "-lfftw3_threads -lfftw3 -lm")

    INCLUDE_DIRECTORIES(/usr/include)
    include_directories(${HDF5_INCLUDE_DIRS})
    include_directories(${CMAKE_SOURCE_DIR}/source/headers)
    find_package(UHD REQUIRED)
    find_package(HDF5 REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS C CXX) # Ensure both C and C++ libraries are found

    include_directories(${UHD_INCLUDE_DIRS})

    LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)
    include_directories(${YOUR_DIRECTORY})
    set(CMAKE_C_OUTPUT_EXTENSION_REPLACE 1)
    set(SOURCES 
                source/src/functions.cpp
                source/src/focus.cpp
    )
    
    add_executable(spectrogram ${SOURCES})
    target_link_libraries(spectrogram ${HDF5_LIBRARIES})

    target_link_libraries(spectrogram HDF5::HDF5)
    target_include_directories(spectrogram PRIVATE ${HDF5_INCLUDE_DIRS})

    set(CPACK_PROJECT_NAME ${PROJECT_NAME})
    set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
    include(CPack)
    find_package(GTest REQUIRED)
    target_link_libraries(spectrogram fftw3)
    include_directories(${GTEST_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(/usr/include)
    include_directories(${HDF5_INCLUDE_DIRS})
    project(spectrogram VERSION 0.1.0 LANGUAGES C CXX)
    target_link_libraries(spectrogram GTest::GTest GTest::Main)
    add_test(NAME ifftTesting COMMAND spectrogram)
    add_test(NAME readHDF COMMAND spectrogram)
    target_link_libraries(spectrogram ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}) # Link HDF5 libraries
endif()