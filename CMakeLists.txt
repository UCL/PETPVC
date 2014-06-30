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

ENABLE_TESTING()
 
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(lib)

ADD_EXECUTABLE(pvc_gtm src/GTM.cxx  )
TARGET_LINK_LIBRARIES(pvc_gtm ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_rbv src/RBV.cxx  )
TARGET_LINK_LIBRARIES(pvc_rbv ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_iy src/IterativeYang.cxx  )
TARGET_LINK_LIBRARIES(pvc_iy ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_mg src/MullerGartner.cxx  )
TARGET_LINK_LIBRARIES(pvc_mg ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_simulate src/Simulate.cxx  )
TARGET_LINK_LIBRARIES(pvc_simulate ${ITK_LIBRARIES})

ADD_SUBDIRECTORY(test)
