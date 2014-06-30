#=========================================================================
#  Author:      Kris Thielemans
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

ADD_EXECUTABLE(pvc_createTestImage CreateTestImage.cxx  )
TARGET_LINK_LIBRARIES(pvc_createTestImage ${ITK_LIBRARIES})

ADD_EXECUTABLE(pvc_compareImages CompareImages.cxx  )
TARGET_LINK_LIBRARIES(pvc_compareImages ${ITK_LIBRARIES})

# There's really only 2 tests currently:
# run IterativeYang and RBV and check that the output is almost 
# identical to the original.
#
# However, to do this, we need to have an input image, and blur it.
# At present, we just do that as seperate tests.
#
# TODO
# KT isn't sure how to make sure that the Compare* tests only run when
# the previous ones were successful.
# This means that they'll fail when using parallel building

ADD_TEST(NAME generate
    COMMAND pvc_createTestImage original.nii mask.nii)

ADD_TEST(NAME simulate
    COMMAND pvc_simulate -x 5 -y 6 -z 7 original.nii filtered.nii )

# macro for adding tests for pvc
# will run matching, underestimated and overestimated blurring
# t_* arguments are thresholds to use for the comparison
macro(add_pvc_test fullname name t_match t_under t_over)
  ADD_TEST(NAME Run${fullname}
    COMMAND pvc_${name} -x 5 -y 6 -z 7 filtered.nii mask.nii ${name}.nii )

  ADD_TEST(NAME Compare${fullname}
    COMMAND pvc_compareImages ${name}.nii original.nii ${t_match} )

  ADD_TEST(NAME Run${fullname}Undercorrect
    COMMAND pvc_${name} -x 4 -y 5 -z 5 filtered.nii mask.nii ${name}_undercorrect.nii )

  ADD_TEST(NAME Compare${fullname}Undercorrect
    COMMAND pvc_compareImages ${name}_undercorrect.nii original.nii ${t_under} )

  ADD_TEST(NAME Run${fullname}Overcorrect
    COMMAND pvc_${name} -x 6 -y 7 -z 8 filtered.nii mask.nii ${name}_overcorrect.nii )

  ADD_TEST(NAME Compare${fullname}Overcorrect
    COMMAND pvc_compareImages ${name}_overcorrect.nii original.nii ${t_over} )
ENDMACRO(add_pvc_test)


add_pvc_test(IterativeYang iy .01 .5 .5 )
add_pvc_test(RBV rbv .01 .5 .5)
