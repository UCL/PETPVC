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

Installation instructions
--------------------------

The following are required to build this software:

	- CMake (http://www.cmake.org/)
	- ITK - Segmentation & Registration Toolkit (http://www.itk.org) >=4.0
	- A C++ compiler

Building for Linux / Unix:

	Create a build directory e.g. :

		mkdir BUILD	

	Change to the build directory :

		cd BUILD

	Run CMake :
	
		ccmake /path/to/src

	Set the build type CMAKE_BUILD_TYPE to 'Release'.

	Generate the project files.

	To build the executables run from the build directory:

		make


Building for Windows:

	Run the CMake GUI :

		cmake-gui

	Specify the source directory as the path to the directory this file is in.

	Choose or create a build directory. 

	Set the build type CMAKE_BUILD_TYPE to 'Release'.

	Generate project files for the appropriate compiler. This project has been 
	tested with Visual Studio 2010 Express. 

	Build the executables accordingly.

Example: Running Iterative Yang PVC
---------------------------

When the project has been successfully built, there will be an executable called 
pvc_iy in the build directory. The program is run using the following command:

	pvc_iy	<petfile> <maskfile> <outputfile> -f <X> <Y> <Z> [-i <V>]

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


In addition, pvc_iy will produce an XML file that conforms to the Slicer execution
model, when the argument "--xml" is specified. The application can therefore be
run from inside 3DSlicer.

This program implements the Iterative Yang PVC technique. 
Please cite the following paper:
	Erlandsson, K. and Buvat, I. and Pretorius, P.H. and 
        Thomas, B.A. and Hutton, B.F., (2012).
	"A review of partial volume correction techniques for emission tomography 
        and their applications in neurology, cardiology and oncology", 
	Physics in Medicine and Biology, vol. 57, no. 21, R119-59.


Contact
--------

Ben Thomas (admin@pet-mri-pvc.org)




