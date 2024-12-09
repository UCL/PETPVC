/*
   petpvcRegionConvolutionPVCImageFilter.h

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

#ifndef __PETPVCREGIONCONVOLUTIONIMAGEFILTER_H
#define __PETPVCREGIONCONVOLUTIONIMAGEFILTER_H

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkImageDuplicator.h>

#include <algorithm>

using namespace itk;

namespace petpvc
{
template< class TInputImage, typename TMaskImage>
class RegionConvolutionPVCImageFilter:public ImageToImageFilter< TInputImage, TInputImage >
{
public:
    /** Standard class typedefs. */
    typedef RegionConvolutionPVCImageFilter             Self;
    typedef ImageToImageFilter< TInputImage, TInputImage > Superclass;
    typedef SmartPointer< Self >        Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
#if ITK_VERSION_MAJOR < 6
    itkTypeMacro(RegionConvolutionPVCImageFilter, ImageToImageFilter);
#else
    itkOverrideGetNameOfClassMacro(RegionConvolutionPVCImageFilter);
#endif

    /** Image related typedefs. */
    typedef TInputImage             InputImageType;
    typedef typename TInputImage::ConstPointer    InputImagePointer;
    typedef typename TInputImage::RegionType RegionType;
    typedef typename TInputImage::SizeType   SizeType;
    typedef typename TInputImage::IndexType  IndexType;
    typedef typename TInputImage::PixelType  PixelType;

    /** Mask image related typedefs. */
    typedef TMaskImage                      MaskImageType;
    typedef typename TMaskImage::ConstPointer    MaskImagePointer;
    typedef typename TMaskImage::RegionType MaskRegionType;
    typedef typename TMaskImage::SizeType   MaskSizeType;
    typedef typename TMaskImage::IndexType  MaskIndexType;
    typedef typename TMaskImage::PixelType  MaskPixelType;

    typedef itk::MultiplyImageFilter<TInputImage, TInputImage> MultiplyFilterType;
    typedef itk::DivideImageFilter<TInputImage,TInputImage, TInputImage> DivideFilterType;
    typedef itk::DiscreteGaussianImageFilter<TInputImage, TInputImage> BlurringFilterType;
    typedef itk::ImageDuplicator<TInputImage> DuplicatorType;

    typedef itk::Vector<float, 3> ITKVectorType;

    /** Image related typedefs. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        3);

    itkStaticConstMacro(MaskImageDimension, unsigned int,
                        3);

    typedef vnl_vector<float> VectorType;
    typedef vnl_matrix<float> MatrixType;

    /** Set the mask image */
    void SetMaskInput(const TMaskImage *input) {
        // Process object is not const-correct so the const casting is required.
        this->SetNthInput( 1, const_cast< TMaskImage * >( input ) );
    }

    /** Get the label image */
    const MaskImageType * GetMaskInput() const {
        return itkDynamicCastInDebugMode< MaskImageType * >( const_cast< DataObject * >( this->ProcessObject::GetInput(0) ) );
    }

    void SetPSF(ITKVectorType vec) {
        this->m_vecVariance = vec;
    }

    ITKVectorType GetPSF() {
        return this->m_vecVariance;
    }

    void SetVerbose( bool bVerbose ) {
        this->m_bVerbose = bVerbose;
    }


protected:
    RegionConvolutionPVCImageFilter();
    ~RegionConvolutionPVCImageFilter() {};

    /** Does the real work. */
    virtual void GenerateData() ITK_OVERRIDE;

    ITKVectorType m_vecVariance;
    bool m_bVerbose;

private:
    RegionConvolutionPVCImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented



};
} //namespace petpvc


#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcRegionConvolutionImageFilter.txx"
#endif


#endif // __PETPVCREGIONCONVOLUTIONIMAGEFILTER_H
