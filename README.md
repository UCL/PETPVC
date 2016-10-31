# PETPVC
PETPVC: toolbox for partial volume correction (PVC) in positron emission tomography (PET)

##Publications
When using this toolbox, please include a reference to the paper:
- <i>PETPVC: a toolbox for performing partial volume correction techniques in positron emission tomography</i><br/>BA Thomas, V Cuplov, A Bousse, A Mendes, K Thielemans, BF Hutton, K Erlandsson<br/>Physics in Medicine and Biology 61 (22), 7975. [DOI](http://dx.doi.org/10.1088/0031-9155/61/22/7975)

---
##Installation instructions

The following are required to build this software:

- [CMake](http://www.cmake.org/)

- [ITK - Segmentation & Registration Toolkit](http://www.itk.org) >=4.7

- A C++ compiler
	
###Building and installing
- Clone the repository
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
##Notes on input and output files
The applications in this toolbox use ITK image readers and writers and can
therefore accept common medical imaging formats such as Nifti, ANALYZE and Nrrd, 
and raw data with an associated meta-data header ([mhd](http://www.itk.org/Wiki/ITK/MetaIO/Documentation#ITK_MetaIO)) file.

The tissue classification maps (referred to as mask files) can either be binary
or probabilistic. All voxel values in a 3D volume must be 0 <= x <= 1. The
PVC applications expect the mask file to be input as a single 4-D volume, where
each 3-D volume consists of a single segmented region. 

The use of 4-D volumes facilitates the use of probabilistic segmentations during 
the PVC. In addition to the constraint that all voxels must be <= 1,  The sum of
a voxel location across the fourth dimension should be <= 1. Ideally it should
be 1, which requires the background to be included as a segmented region.

###Special cases where the inputs/outputs are different
####Muller-Gartner (MG):
The Muller-Gartner correction requires only the grey matter and white
matter masks. Technically, the CSF space should be included as a third 
region, but the contribution of this region is assumed to be zero. The 
MG application still requires a 4-D mask volume, where the first volume is
grey matter and the second is white matter. The order is important. The 
4-D mask file can contain more than two 3-D volumes, but these will be 
ignored by the MG PVC.

####Geometric Transfer Matrix (GTM) method:
GTM cannot produce an image. The output of the GTM is a 
comma-separated value (CSV) file of regional mean values. The order of the
mean values for each region is written in the same order as they appear in
the fourth dimension of the mask file.
