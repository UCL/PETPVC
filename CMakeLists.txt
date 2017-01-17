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

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)

MAKE_DIRECTORY(${PROJECT_BINARY_DIR}/config)

CONFIGURE_FILE(
        "${PROJECT_SOURCE_DIR}/EnvironmentInfo.h.in"
        "${PROJECT_BINARY_DIR}/config/EnvironmentInfo.h"
)

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/config")

ENABLE_TESTING()
 
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(lib)

ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(src)

SET(CPACK_PACKAGE_NAME "PETPVC")
SET(CPACK_PACKAGE_VENDOR "Benjamin A. Thomas")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Installer for PETPVC")
SET(CPACK_PACKAGE_VERSION_MAJOR VERSION_MAJOR)
SET(CPACK_PACKAGE_VERSION_MINOR VERSION_MINOR)
SET(CPACK_PACKAGE_VERSION_PATCH VERSION_PACK)
SET(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Installer")

SET(CPACK_COMPONENTS_ALL apps)
SET(CPACK_COMPONENT_APPS_DISPLAY_NAME "PETPVC applications")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")

INCLUDE(CPack)