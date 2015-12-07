/*
   petpvcRLPVCImageFilter.h

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

#ifndef __PETPVCRLPVCIMAGEFILTER_H
#define __PETPVCRLPVCIMAGEFILTER_H

#include "itkImage.h"
#include "itkImageToImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkImageDuplicator.h>

#include <algorithm>

using namespace itk;

namespace petpvc
{
template< class TInputImage>
class RichardsonLucyPVCImageFilter:public ImageToImageFilter< TInputImage, TInputImage >
{
public:
    /** Standard class typedefs. */
    typedef RichardsonLucyPVCImageFilter             Self;
    typedef ImageToImageFilter< TInputImage, TInputImage > Superclass;
    typedef SmartPointer< Self >        Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(RichardsonLucyPVCImageFilter, ImageToImageFilter);

    /** Image related typedefs. */
    typedef TInputImage             InputImageType;
    typedef typename TInputImage::ConstPointer    InputImagePointer;
    typedef typename TInputImage::RegionType RegionType;
    typedef typename TInputImage::SizeType   SizeType;
    typedef typename TInputImage::IndexType  IndexType;
    typedef typename TInputImage::PixelType  PixelType;

	typedef itk::ImageRegionConstIterator<TInputImage> ConstImageIterator;

    //For calculating mean values from image
    typedef itk::StatisticsImageFilter<TInputImage> StatisticsFilterType;
    //Extracts a 3D volume from 4D file.
 
    typedef itk::MultiplyImageFilter<TInputImage, TInputImage> MultiplyFilterType;
	typedef itk::DivideImageFilter<TInputImage, TInputImage, TInputImage> DivideFilterType;
    typedef itk::AddImageFilter<TInputImage, TInputImage> AddFilterType;
	typedef itk::SubtractImageFilter<TInputImage, TInputImage> SubFilterType;
    typedef itk::DiscreteGaussianImageFilter<TInputImage, TInputImage> BlurringFilterType;
	typedef itk::ThresholdImageFilter<TInputImage> ThresholdFilterType;
    typedef itk::ImageDuplicator<TInputImage> DuplicatorType;

    typedef itk::Vector<float, 3> ITKVectorType;

    /** Image related typedefs. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        3);

    void SetPSF(ITKVectorType vec) {
        this->m_vecVariance = vec;
    }


    ITKVectorType GetPSF() {
        return this->m_vecVariance;
    }

    void SetIterations( unsigned int nIters ) {
        this->m_nIterations = nIters;
    }

    void SetStoppingCond( float stop ) {
        this->m_fStopCriterion = stop;
    }

    void SetVerbose( bool bVerbose ) {
        this->m_bVerbose = bVerbose;
    }


protected:
    RichardsonLucyPVCImageFilter();
    ~RichardsonLucyPVCImageFilter() {};

    /** Does the real work. */
    virtual void GenerateData();

    ITKVectorType m_vecVariance;
    unsigned int m_nIterations;
    float m_fStopCriterion;
    bool m_bVerbose;

private:
    RichardsonLucyPVCImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented



};
} //namespace petpvc


#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcRLPVCImageFilter.txx"
#endif


#endif // __PETPVCRichardsonLucyIMAGEFILTER_H
