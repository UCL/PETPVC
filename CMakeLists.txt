#=========================================================================
#  Authors:      Benjamin A. Thomas, Kris Thielemans
#
#  Copyright University College London
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#=========================================================================

cmake_minimum_required(VERSION 2.8)

PROJECT(PETPVC)

# set default build-type to Release
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 2)
SET(VERSION_PATCH 10)

MAKE_DIRECTORY(${PROJECT_BINARY_DIR}/config)

CONFIGURE_FILE(
        "${PROJECT_SOURCE_DIR}/EnvironmentInfo.h.in"
        "${PROJECT_BINARY_DIR}/config/EnvironmentInfo.h"
)

#if(CMAKE_COMPILER_IS_GNUCXX)
#    include(CodeCoverage)
#    setup_target_for_coverage(${PROJECT_NAME}_coverage tests coverage)
#endif()

#SET(CMAKE_CXX_FLAGS "-Wall -pedantic -pthread -g -O0 --coverage")

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/config")

ENABLE_TESTING()

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

IF(NOT ITKReview_LOADED)
        MESSAGE(FATAL_ERROR "ITK must be built with Module_ITKReview enabled.")
ENDIF()

INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(lib)

ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(src)

INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/parc/" DESTINATION "parc")

SET(CPACK_PACKAGE_NAME "PETPVC")
SET(CPACK_PACKAGE_VENDOR "Institute of Nuclear Medicine")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Installer for PETPVC")
SET(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
SET(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "PETPVC")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")

SET(CPACK_COMPONENTS_ALL applications)
SET(CPACK_COMPONENT_APPS_DISPLAY_NAME "PETPVC applications")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")

INCLUDE(CPack)
