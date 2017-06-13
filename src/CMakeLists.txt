#=========================================================================
#  Author:      Benjamin A. Thomas
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

ADD_EXECUTABLE(pvc_gtm GTM.cxx  )
TARGET_LINK_LIBRARIES(pvc_gtm ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_labbe Labbe.cxx  )
TARGET_LINK_LIBRARIES(pvc_labbe ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_mtc MTC.cxx  )
TARGET_LINK_LIBRARIES(pvc_mtc ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_rbv RBV.cxx  )
TARGET_LINK_LIBRARIES(pvc_rbv ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_iy IterativeYang.cxx  )
TARGET_LINK_LIBRARIES(pvc_iy ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_diy DiscreteIY.cxx  )
TARGET_LINK_LIBRARIES(pvc_diy ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_mg MullerGartner.cxx  )
TARGET_LINK_LIBRARIES(pvc_mg ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_vc VanCittert.cxx  )
TARGET_LINK_LIBRARIES(pvc_vc ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_rl RichardsonLucy.cxx  )
TARGET_LINK_LIBRARIES(pvc_rl ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_stc STC.cxx  )
TARGET_LINK_LIBRARIES(pvc_stc ${ITK_LIBRARIES})

ADD_EXECUTABLE(petpvc PETPVC.cxx  )
TARGET_LINK_LIBRARIES(petpvc ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_simulate Simulate.cxx  )
TARGET_LINK_LIBRARIES(pvc_simulate ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_make4d Make4D.cxx  )
TARGET_LINK_LIBRARIES(pvc_make4d ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_relabel Relabel.cxx  )
TARGET_LINK_LIBRARIES(pvc_relabel ${ITK_LIBRARIES})

INSTALL(TARGETS pvc_gtm pvc_labbe pvc_rbv pvc_iy
				pvc_diy pvc_mg pvc_mtc pvc_vc pvc_rl pvc_stc 
				petpvc pvc_simulate pvc_make4d
				pvc_relabel DESTINATION bin
				COMPONENT applications)
