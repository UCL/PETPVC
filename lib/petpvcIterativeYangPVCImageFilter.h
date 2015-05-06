/*
   petpvcIterativeYangPVCImageFilter.h

   Author:      Benjamin A. Thomas

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

#ifndef __PETPVCITERATIVEYANGIMAGEFILTER_H
#define __PETPVCITERATIVEYANGIMAGEFILTER_H

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "petpvcFuzzyCorrectionFilter.h"

#include <itkExtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkImageDuplicator.h>

#include <algorithm>

using namespace itk;

namespace petpvc
{
template< class TInputImage, typename TMaskImage>
class IterativeYangPVCImageFilter:public ImageToImageFilter< TInputImage, TInputImage >
{
public:
    /** Standard class typedefs. */
    typedef IterativeYangPVCImageFilter             Self;
    typedef ImageToImageFilter< TInputImage, TInputImage > Superclass;
    typedef SmartPointer< Self >        Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(IterativeYangPVCImageFilter, ImageToImageFilter);

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

    //For calculating mean values from image
    typedef itk::StatisticsImageFilter<TInputImage> StatisticsFilterType;
    //Extracts a 3D volume from 4D file.
    typedef itk::ExtractImageFilter<TMaskImage, TInputImage> ExtractFilterType;
    typedef itk::MultiplyImageFilter<TInputImage, TInputImage> MultiplyFilterType;
    typedef itk::DivideImageFilter<TInputImage,TInputImage, TInputImage> DivideFilterType;
    typedef itk::AddImageFilter<TInputImage, TInputImage> AddFilterType;
    typedef itk::DiscreteGaussianImageFilter<TInputImage, TInputImage> BlurringFilterType;
    typedef itk::ImageDuplicator<TInputImage> DuplicatorType;

    typedef FuzzyCorrectionFilter<TMaskImage> FuzzyCorrFilterType;
    typedef itk::Vector<float, 3> ITKVectorType;

    /** Image related typedefs. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        3);

    itkStaticConstMacro(MaskImageDimension, unsigned int,
                        4);

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

    VectorType GetCorrectedMeans() const {
        return this->m_vecRegMeansPVCorr;
    }

    MatrixType GetMatrix() const {
        return this->m_matGTM;
    }

    void SetPSF(ITKVectorType vec) {
        this->m_vecVariance = vec;
    }


    ITKVectorType GetPSF() {
        return this->m_vecVariance;
    }

    void SetIterations( unsigned int nIters ) {
        this->m_nIterations = nIters;
    }

    void SetVerbose( bool bVerbose ) {
        this->m_bVerbose = bVerbose;
    }


protected:
    IterativeYangPVCImageFilter();
    ~IterativeYangPVCImageFilter() {};

    /** Does the real work. */
    virtual void GenerateData();

    VectorType m_vecRegMeansPVCorr;
    MatrixType m_matGTM;
    ITKVectorType m_vecVariance;
    unsigned int m_nIterations;
    bool m_bVerbose;

private:
    IterativeYangPVCImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented



};
} //namespace petpvc


#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcIterativeYangPVCImageFilter.txx"
#endif


#endif // __PETPVCITERATIVEYANGIMAGEFILTER_H
