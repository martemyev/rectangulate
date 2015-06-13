cmake_minimum_required(VERSION 2.8)


#-------------------------------------------------------------------------------
# Build type.
#
# There are three of them: Release (default), Debug, and RelWithDebInfo. See the
# section 'Choosing the build type' below for details.
#-------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build: Debug | Release | RelWithDebInfo" FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)


#-------------------------------------------------------------------------------
# Project.
#
# Setting the project name and some auxiliary compilation flags
#-------------------------------------------------------------------------------
project(rectangulate)
message(STATUS "C++ compiler ID = ${CMAKE_CXX_COMPILER_ID}")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(MY_CXX_FLAGS "-std=c++0x -Wall -Wextra -Wshadow -Weffc++ -pedantic -Werror")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(MY_CXX_FLAGS "-std=c++0x -Wall -Wextra -Wshadow -Weffc++ -pedantic -Werror")
else() 
  set(MY_CXX_FLAGS "-std=c++0x -Wall -pedantic -Werror")
endif()


#-------------------------------------------------------------------------------
# Version of the project.
#-------------------------------------------------------------------------------
set(VERSION_MAJOR 0)
set(VERSION_MINOR 0)
set(VERSION_PATCH 1)


#-------------------------------------------------------------------------------
# Choosing the build type.
#
# The default build type is RELEASE. It is good for running the simulation,
# but very limited in checking errors while executing. Therefore once the
# program stops unexpectedly or produces strange results, one should first
# of all rebuild this this project in the DEBUG mode. In this case much more
# errors checks are performed revealing errors, if any. Be aware, however,
# that the program runs much slower in this mode.
# There is also a RelWithDebInfo mode which stands for "release mode with
# debugging information". In this mode the project is built with both an
# optimization key '-O2' and with debugging symbols '-g'. This is ideal for
# using in some of the Valgrind programs (for example, callgrind). However, if
# one intends to check memory leaks with the memcheck program, they should use
# the DEBUG version of the project. We also add -fno-inline option in non
# release modes for better Valgrind profiling.
#-------------------------------------------------------------------------------
if(BUILD_TYPE MATCHES DEBUG)
  set(DEBUG_REC ON CACHE BOOL
      "This option is for internal usage - not the compiler debug options" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${MY_CXX_FLAGS} -fno-inline")

elseif(BUILD_TYPE MATCHES RELWITHDEBINFO)
  set(DEBUG_REC OFF CACHE BOOL
      "This option is for internal usage - not the compiler debug options" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${MY_CXX_FLAGS} -fno-inline")

elseif(BUILD_TYPE MATCHES RELEASE)
  set(DEBUG_REC OFF CACHE BOOL
      "This option is for internal usage - not the compiler debug options" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MY_CXX_FLAGS}")

else()
  message(FATAL_ERROR "Unknown build type: ${BUILD_TYPE}. Options: Debug | Release | RelWithDebInfo")
endif()

set(SELFCHECK OFF CACHE BOOL
    "Functions check themselves - it slows the code down, but it's useful for debugging" FORCE)


#-------------------------------------------------------------------------------
# 64-bitness
#-------------------------------------------------------------------------------
set(HAVE_64BIT_SIZE_T OFF CACHE INTERNAL "")
include(CheckTypeSize)
check_type_size("void*" SIZEOF_VOID_P)
if(SIZEOF_VOID_P EQUAL 8)
  set(HAVE_64BIT_SIZE_T ON CACHE INTERNAL "" FORCE)
endif(SIZEOF_VOID_P EQUAL 8)


#-------------------------------------------------------------------------------
# Config file
#-------------------------------------------------------------------------------
configure_file(
  "${PROJECT_SOURCE_DIR}/config.hpp.in"
  "${PROJECT_SOURCE_DIR}/headers/config.hpp")


#-------------------------------------------------------------------------------
# Headers and sources
#-------------------------------------------------------------------------------
FILE(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/sources/*.cpp")
FILE(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/headers/*.hpp")
include_directories(${PROJECT_SOURCE_DIR}/headers)


#-------------------------------------------------------------------------------
# Put the information about the configuration on the screen
#-------------------------------------------------------------------------------
  message("\n--------------- configuration information ----------------")
  message("build type     " ${BUILD_TYPE})
  message("compiler       " ${CMAKE_CXX_COMPILER})
if(BUILD_TYPE MATCHES DEBUG)
  message("compiler flags " ${CMAKE_CXX_FLAGS_DEBUG})
elseif(BUILD_TYPE MATCHES RELWITHDEBINFO)
  message("compiler flags " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
elseif(BUILD_TYPE MATCHES RELEASE)
  message("compiler flags " ${CMAKE_CXX_FLAGS_RELEASE})
endif()

  message("----------------------------------------------------------")
  message("DEBUG_REC         " ${DEBUG_REC})
  message("HAVE_64BIT_SIZE_T " ${HAVE_64BIT_SIZE_T})
  message("----------------------------------------------------------")


#-------------------------------------------------------------------------------
# Declaring the executable
#-------------------------------------------------------------------------------
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})


#-------------------------------------------------------------------------------
# Packaging of the project
#-------------------------------------------------------------------------------
set(PACKAGE_FILE_NAME
    "${PROJECT_NAME}-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

add_custom_target(dist
                  COMMAND git archive --prefix=${PACKAGE_FILE_NAME}/ -o ${PROJECT_BINARY_DIR}/${PACKAGE_FILE_NAME}.tar.gz HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})


#-------------------------------------------------------------------------------
# Further instructions
#-------------------------------------------------------------------------------
message("\n======================= building ==============================")
message("To build the project, launch:")
message("  $ make")
message("===============================================================\n")

