/*
   petpvcVanCittertPVCImageFilter.h

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

#ifndef __PETPVCVANCITTERTIMAGEFILTER_H
#define __PETPVCVANCITTERTIMAGEFILTER_H

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
class VanCittertPVCImageFilter:public ImageToImageFilter< TInputImage, TInputImage >
{
public:
    /** Standard class typedefs. */
    typedef VanCittertPVCImageFilter             Self;
    typedef ImageToImageFilter< TInputImage, TInputImage > Superclass;
    typedef SmartPointer< Self >        Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(VanCittertPVCImageFilter, ImageToImageFilter);

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

    void SetAlpha( float alpha ) {
        this->m_fAlpha = alpha;
    }

    void SetStoppingCond( float stop ) {
        this->m_fStopCriterion = stop;
    }

    void SetVerbose( bool bVerbose ) {
        this->m_bVerbose = bVerbose;
    }


protected:
    VanCittertPVCImageFilter();
    ~VanCittertPVCImageFilter() {};

    /** Does the real work. */
    virtual void GenerateData();

    ITKVectorType m_vecVariance;
    unsigned int m_nIterations;
    float m_fAlpha;
    float m_fStopCriterion;
    bool m_bVerbose;

private:
    VanCittertPVCImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented



};
} //namespace petpvc


#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcVanCittertPVCImageFilter.txx"
#endif


#endif // __PETPVCVANCITTERTIMAGEFILTER_H
