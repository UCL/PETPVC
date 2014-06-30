/*
   PET Partial Volume Correction (PVC) toolbox

   Copyright 2013 Institute of Nuclear Medicine, University College London.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

--------------------------
Installation instructions
--------------------------

The following are required to build this software:

	- CMake (http://www.cmake.org/)
	- ITK - Segmentation & Registration Toolkit (http://www.itk.org) >=4.0
	- A C++ compiler


Building for Linux / Unix /  etc:

step 1: Run CMake

   a) Using the CMake GUI :

		cmake-gui

	Set the build type CMAKE_BUILD_TYPE to 'Release' and choose an installation directory.

	Generate project files for the appropriate compiler. 


   b) Using a terminal-driven CMake (Unix / Linux)

	Create a build directory e.g. :

		mkdir BUILD	

	Change to the build directory :

		cd BUILD

	Run CMake :
	
		ccmake /path/to/src

	Set the build type CMAKE_BUILD_TYPE to 'Release' and choose an installation directory.

	Generate the project files.

step 2: building

	This step depends on which type of "generator" for the project you used.

	If you used "unix makefiles" (default on Linux), you need to build
	the executables run from the build directory:

		make
		make test
		make install

	(If you used a system path for the installation directory, you might need to
	use "sudo" or change user for the installation step.)

	If you created project files for an IDE (e.g. Visual Studio, XCode, Eclipse),
	you will need to open them in your IDE (e.g. PETPVC.sln for Visual Studio) and
        build the targets BUILD_ALL, RUN_TESTS, INSTALL	

The executables will be installed in the bin subdirectory of your chosen installation-prefix.

This project has been tested with gcc, Visual Studio 2010 Express, mingw-w32 and mingw-64.

--------------------------------
Notes on input and output files
--------------------------------

The applications in this toolbox use ITK image readers and writers and can
therefore accept common medical imaging formats such as Nifti, ANALYZE and Nrrd, 
and raw data with an associated meta-data header (mhd) file
(see : http://www.itk.org/Wiki/ITK/MetaIO/Documentation#ITK_MetaIO).

The tissue classification maps (referred to as mask files) can either be binary
or probabilistic. All voxel values in a 3D volume must be 0 <= x <= 1. The
PVC applications expect the mask file to be input as a single 4-D volume, where
each 3-D volume consists of a single segmented region. 

The use of 4-D volumes facilitates the use of probabilistic segmentations during 
the PVC. In addition to the constraint that all voxels must be <= 1,  The sum of
a voxel location across the fourth dimension should be <= 1. Ideally it should
be 1, which requires the background to be included as a segmented region.

Special cases where the inputs/outputs are different
-----------------------------------------------------

    Muller-Gartner (MG):

    The Muller-Gartner correction requires only the grey matter and white
    matter masks. Technically, the CSF space should be included as a third 
    region, but the contribution of this region is assumed to be zero. The 
    MG application still requires a 4-D mask volume, where the first volume is
    grey matter and the second is white matter. The order is important. The 
    4-D mask file can contain more than two 3-D volumes, but these will be 
    ignored by the MG PVC.

    Geometric Transfer Matrix (GTM) method:

    By design the GTM cannot produce an image. The output of the GTM is a 
    comma-separated value (CSV) file of regional mean values. The order of the
    mean values for each region is written in the same order as they appear in
    the fourth dimension of the mask file.

-------------------------    
Application descriptions
-------------------------

Running Iterative Yang (IY) PVC
--------------------------------

When the project has been successfully built, there will be an executable called 
pvc_iy in the build directory. The program is run using the following command:

	pvc_iy	<petfile> <maskfile> <outputfile> -x <X> -y <Y> -z <Z> [-i <V>]

where:

	Compulsory:

	<petfile> is the reconstructed PET emission file.
	<maskfile> is the filename of the segmentation. This a 4D file containing 
		1 segmented region per 3D volume.
	<outputfile> is the resultant PV-corrected PET file to be created.
	<X,Y,Z> are the measure of the full-width at half maximum (FWHM) of the
		scanner PSF in the x-, y- and z-axis, specified in mm. The PSF is
		assumed to be Gaussian, although can be changed in the source 
		code.

	Optional:

	-i <V> controls the number of iterations to be performed. If not supplied, 
		the default value of 5 is used.

In addition, pvc_iy (and the other applications) will produce an XML output that
conforms to the Slicer execution model, when the argument "--xml" is specified.
The application can therefore be run from inside 3DSlicer.

This program implements the Iterative Yang PVC technique. 
The method is described in:
	Erlandsson, K. and Buvat, I. and Pretorius, P.H. and 
        Thomas, B.A. and Hutton, B.F., (2012).
	"A review of partial volume correction techniques for emission tomography 
        and their applications in neurology, cardiology and oncology", 
	Physics in Medicine and Biology, vol. 57, no. 21, R119-59.


Running region-based voxel-wise (RBV) correction
-------------------------------------------------

The executable pvc_rbv will perform the RBV PVC technique. The program is run 
using the following command:

	pvc_rbv	<petfile> <maskfile> <outputfile> -x <X> -y <Y> -z <Z> 

where:

	Compulsory:

	<petfile> is the reconstructed PET emission file.
	<maskfile> is the filename of the segmentation. This a 4D file containing 
		1 segmented region per 3D volume.
	<outputfile> is the resultant PV-corrected PET file to be created.
	<X,Y,Z> are the measure of the full-width at half maximum (FWHM) of the
		scanner PSF in the x-, y- and z-axis, specified in mm. The PSF is
		assumed to be Gaussian, although can be changed in the source 
		code.

This program implements the region-based voxel-wise (RBV) PVC technique. 
The method is described in:

        Thomas, B. and Erlandsson, K. and Modat, M. and Thurfjell, L. and
        Vandenberghe, R. and Ourselin, S. and Hutton, B. (2011). "The importance
        of appropriate partial volume correction for PET quantification in 
        Alzheimer's disease". European Journal of Nuclear Medicine and 
        Molecular Imaging, 38:1104-1119.


Running Geometric Transfer Matrix (GTM) correction
---------------------------------------------------

The executable pvc_gtm will perform the GTM PVC technique. The program is run 
using the following command:

	pvc_gtm	<petfile> <maskfile> <outputfile> -x <X> -y <Y> -z <Z> 

where:

	Compulsory:

	<petfile> is the reconstructed PET emission file.
	<maskfile> is the filename of the segmentation. This a 4D file containing 
		1 segmented region per 3D volume.
	<outputfile> is a comma-separated value (CSV) file containing regional
                mean values. The values are output in the same order as the 
                region masks are defined in the <maskfile>.
	<X,Y,Z> are the measure of the full-width at half maximum (FWHM) of the
		scanner PSF in the x-, y- and z-axis, specified in mm. The PSF is
		assumed to be Gaussian, although can be changed in the source 
		code.

This program implements the Geometric Transfer Matrix (GTM) PVC technique. 
The method is described in:
        Rousset, O. G. and Ma, Y. and Evans, A. C. (1998). "Correction for 
        partial volume effects in PET: principle and validation". Journal of 
        Nuclear Medicine, 39(5):904-11.


Running Muller-Gartner (MG) correction
---------------------------------------

The executable pvc_mg will perform the MG PVC technique. The program is run 
using the following command:

	pvc_mg	<petfile> <maskfile> <outputfile> -x <X> -y <Y> -z <Z> 

where:

	Compulsory:

	<petfile> is the reconstructed PET emission file.
	<maskfile> is the filename of the segmentation. This a 4D file containing 
		1 segmented region per 3D volume. The first region must be the
                the grey matter mask, the second must be the white matter mask.
	<outputfile> is a comma-separated value (CSV) file containing regional
                mean values. The values are output in the same order as the 
                region masks are defined in the <maskfile>.
	<X,Y,Z> are the measure of the full-width at half maximum (FWHM) of the
		scanner PSF in the x-, y- and z-axis, specified in mm. The PSF is
		assumed to be Gaussian, although can be changed in the source 
		code.

This program implements the Muller-Gartner (MG) PVC technique. 
The method is described in:
        Muller-Gartner, H. W. et al. (1992). "Measurement of radiotracer 
        concentration in brain gray matter using positron emission 
        tomography: MRI-based correction for partial volume effects". 
        J Cereb Blood Flow Metab, 12(4), 571-83.

--------
Contact
--------

Ben Thomas (admin@pet-mri-pvc.org)




