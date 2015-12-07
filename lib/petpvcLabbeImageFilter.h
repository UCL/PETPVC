/*
   petpvcLabbeImageFilter.h

   Author:      Benjamin A. Thomas

   Copyright 2015 Institute of Nuclear Medicine, University College London.
   Copyright 2015 Clinical Imaging Research Centre, A*STAR-NUS.

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

#ifndef __PETPVCLABBEIMAGEFILTER_H
#define __PETPVCLABBEIMAGEFILTER_H

#include "itkImageToImageFilter.h"
#include "vnl/vnl_matrix.h"
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <itkImage.h>

//A class to perform Labbe PVC.

using namespace itk;


namespace petpvc
{

template<class TImage>
class LabbeImageFilter : public ImageToImageFilter<TImage, TImage>
{
public:

    typedef LabbeImageFilter Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self> Pointer;

    //Matrix to hold correction factors
    typedef vnl_matrix<float> MatrixType;

    //Vector containing size of region.
    typedef vnl_vector<float> VectorType;
    typedef itk::Vector<float, 3> ITKVectorType;

    itkNewMacro(Self);

    itkTypeMacro(LabbeImageFilter, ImageToImageFilter);

    //Returns correction factors.
    vnl_matrix<float> GetMatrix() {
        return *this->matCorrFactors;
    };

    //Returns region size.

    vnl_vector<float> GetSumOfRegions() {
        return *this->vecSumOfRegions;
    };

    void SetPSF(ITKVectorType vec) {
        this->vecVariance = vec;
    };

    ITKVectorType GetPSF() {
        return this->vecVariance;
    };


protected:
    LabbeImageFilter();
    ~LabbeImageFilter();

    /** Does the real work. */
    virtual void GenerateData();

private:
    LabbeImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &); //purposely not implemented

    MatrixType *matCorrFactors;
    VectorType *vecSumOfRegions;
    ITKVectorType vecVariance;
};
} //namespace PETPVC

#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcLabbeImageFilter.txx"
#endif

#endif // __PETPVCLabbeIMAGEFILTER_H
