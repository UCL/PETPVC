# PETPVC 
[![Build and ctest and recon_test_pack CI](https://github.com/UCL/PETPVC/actions/workflows/build-test.yml/badge.svg)](https://github.com/UCL/PETPVC/actions/workflows/build-test.yml)[![Build Status](https://dev.azure.com/petpvc/petpvc/_apis/build/status/UCL.PETPVC?branchName=master)](https://dev.azure.com/petpvc/petpvc/_build/latest?definitionId=1&branchName=master) [![Build status](https://ci.appveyor.com/api/projects/status/7kk9ua9r0lybinwa/branch/master?svg=true)](https://ci.appveyor.com/project/bathomas/petpvc/branch/master) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/ab83c41f20194c2c82fbc74e8788f928)](https://www.codacy.com/app/bathomas/PETPVC?utm_source=github.com&utm_medium=referral&utm_content=UCL/PETPVC&utm_campaign=badger) [![DOI](https://zenodo.org/badge/17082200.svg)](https://zenodo.org/badge/latestdoi/17082200)

PETPVC: toolbox for partial volume correction (PVC) in positron emission tomography (PET)

## Publications
When using this toolbox, please include a reference to the paper:
- <i>PETPVC: a toolbox for performing partial volume correction techniques in positron emission tomography</i><br/>BA Thomas, V Cuplov, A Bousse, A Mendes, K Thielemans, BF Hutton, K Erlandsson<br/>Physics in Medicine and Biology 61 (22), 7975. [DOI](http://dx.doi.org/10.1088/0031-9155/61/22/7975)

---
## Pre-built binaries

Binaries for Linux, Mac and Windows are provided in the ['Releases' section on Github](https://github.com/UCL/PETPVC/releases). When running the Windows version, you must have installed the ***Visual C++ Redistributable Packages for Visual Studio 2019*** ([link](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)). For 64-bit Windows, you might have to install both ```vcredist_x86.exe``` and ```vcredist_x64.exe```.

Alternatively, PETPVC can also be installed via `conda`, see https://anaconda.org/conda-forge/petpvc.

---
## Installation from source instructions

The following are required to build this software:

- [CMake](http://www.cmake.org/)

- [ITK - Segmentation & Registration Toolkit](http://www.itk.org) >=4.7

- A C++ compiler
	
### Building and installing
- Ensure that ```ITK``` has been built successfully, with the ```ITKReview``` module (```Module_ITKReview``` in CMake) enabled.
- Clone this repository
```bash
git clone https://github.com/UCL/PETPVC.git
```
- Create a build directory
```bash
mkdir BUILD
```
- Change to the build directory
```bash
cd BUILD
```
- Run CMake
```bash
cmake /path/to/repository
```
- Build and install
```bash
make
make test
make install
```
---

## Usage

### PVC methods

Although there are currently executables for every method, it is based to use
`petpvc` which has options to specify the method and parameters. Type
`petpvc` without arguments to get a usage message.

An example of running iterative Yang with a 6mm PSF:

```
	petpvc -i <PET> -m <MASK> -o <OUTPUT> --pvc IY -x 6.0 -y 6.0 -z 6.0 [--debug]
```
where ```<PET>``` is the PET image file, ```<MASK>``` is the 4-D mask image file and ```<OUTPUT>``` is the destination file for the PV-corrected image.

*Warning*: there are currently 2 options which seem the same:
- `-n`: specifies the number of iterations for iterative Yang
- `-k`: specifies the number of iterations for deconvolution methods such as
Van Cittert and Richardson-Lucy.

Therefore, if you use RL only, you have to use the `-k` option. You can use
both options with for instance `IY+RL` to first run iterative Yang followed by
Richardson-Lucy for extra deconvolution.

### Extras

In addition, there are some utilities that you might find useful:
- `pvc_simulate` allows you to blur an image with a Gaussian (e.g. to simulate
resolution effects)
- some [mask related tools](parc/README.md)

---
## Notes on input and output files
The applications in this toolbox use ITK image readers and writers and can
therefore accept common medical imaging formats such as Nifti, ANALYZE and Nrrd, and raw data with an associated meta-data header ([mhd](http://www.itk.org/Wiki/ITK/MetaIO/Documentation#ITK_MetaIO)) file.

The tissue classification maps (referred to as mask files) can either be binary or probabilistic. All voxel values in a 3-D volume must be 0 <= x <= 1. The PVC applications expect the mask file to be input as a single 4-D volume, where each 3-D volume consists of a single segmented region. 

The use of 4-D volumes facilitates the use of probabilistic segmentations during the PVC. In addition to the constraint that all voxels must be <= 1,  The sum of a voxel location across the fourth dimension should be <= 1. Ideally it should be 1, which requires the background to be included as a segmented region.

### Special cases where the inputs/outputs are different
#### Muller-Gartner (MG):
The Muller-Gartner correction requires only the grey matter and white
matter masks. Technically, the CSF space should be included as a third 
region, but the contribution of this region is assumed to be zero. The 
MG application still requires a 4-D mask volume, where the first volume is
grey matter and the second is white matter. The order is important. The 
4-D mask file can contain more than two 3-D volumes, but these will be 
ignored by the MG PVC.

#### Geometric Transfer Matrix (GTM) method:
GTM cannot produce an image. The output of the GTM is a 
comma-separated value (CSV) file of regional mean values. The order of the
mean values for each region is written in the same order as they appear in
the fourth dimension of the mask file.

#### Single Target Correction (STC) method:
The STC method corrects a single region. The mask image should be a 3-D volume, where each voxel in the target region should be 1. All other voxels should be 0. 
