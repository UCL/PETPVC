/*
   petpvcFuzzyCorrectionFilter.h

   Author:      Benjamin A. Thomas

   Copyright 2013-2015 Institute of Nuclear Medicine, University College London.
   Copyright 2014-2015 Clinical Imaging Research Centre, A*STAR-NUS.

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

#ifndef __PETPVCFUZZYCORRECTIONFILTER_H
#define __PETPVCFUZZYCORRECTIONFILTER_H

#include "itkImageToImageFilter.h"
#include "vnl/vnl_matrix.h"

//A class to perform fuzziness 'correction'. This is required to weight the
//mean values correctly when using probabilistic segmentations. For binary
//(piece-wise constant) segmentations, this is not necessary and will return
//the identity matrix.

using namespace itk;

namespace petpvc
{

template<class TImage>
class FuzzyCorrectionFilter : public ImageToImageFilter<TImage, TImage>
{
public:

    typedef FuzzyCorrectionFilter Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self> Pointer;

    //Matrix to hold correction factors
    typedef vnl_matrix<float> MatrixType;

    //Vector containing size of region.
    typedef vnl_vector<float> VectorType;

    itkNewMacro(Self);

    itkTypeMacro(FuzzyCorrFilter, ImageToImageFilter);

    //Returns correction factors.

    vnl_matrix<float> GetMatrix() {
        return *this->matFuzz;
    };

    //Returns region size.

    vnl_vector<float> GetSumOfRegions() {
        return *this->vecSumOfRegions;
    };

protected:
    FuzzyCorrectionFilter();
    ~FuzzyCorrectionFilter();

    /** Does the real work. */
    virtual void GenerateData();


private:
    FuzzyCorrectionFilter(const Self &); //purposely not implemented
    void operator=(const Self &); //purposely not implemented

    MatrixType *matFuzz;
    VectorType *vecSumOfRegions;

};
} //namespace PETPVC

#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcFuzzyCorrectionFilter.txx"
#endif

#endif // __PETPVCFUZZYCORRECTIONFILTER_H
