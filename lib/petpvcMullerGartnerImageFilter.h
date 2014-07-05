/*
   petpvcMullerGartnerImageFilter.h

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

   This program implements the Muller-Gartner (MG) partial volume
   correction (PVC) technique. The method is described in:
        Muller-Gartner, H. W. et al. (1992). Measurement of radiotracer
                concentration in brain gray matter using positron emission
                tomography: MRI-based correction for partial volume effects.
                J Cereb Blood Flow Metab, 12(4), 571-83.

 */
#ifndef __PETPVCMULLERGARTNERIMAGEFILTER_H
#define __PETPVCMULLERGARTNERIMAGEFILTER_H

#include "itkImage.h"
#include "itkInPlaceImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

using namespace itk;

#ifndef VectorType
typedef itk::Vector<float, 3> VectorType;
#endif

namespace petpvc
{

/** \class MullerGartnerImageFilter
 *
 * \brief Performs Muller-Gartner partial-volume correction.
 *
 * This filter performs Muller-Gartner partial-volume correction for PET.
 * Registered and segmented MR brain tissue images are required in addition
 * to the PET. An accurate measurement of the point-spread function (PSF)
 * in mm is necessary for correction.
 *
 * \author Benjamin A. Thomas, Institute of Nuclear Medicine, UCLH (2009)
 *
 */

template <class TInputImage1, class TInputImage2, class TInputImage3,
         class TOutputImage>
class MullerGartnerImageFilter :
    public InPlaceImageFilter<TInputImage1, TOutputImage>
{
public:

    typedef MullerGartnerImageFilter Self;
    typedef InPlaceImageFilter<TInputImage1, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self);

    itkTypeMacro(MullerGartnerImageFilter, InPlaceImageFilter);

    //PET image
    typedef TInputImage1 Input1ImageType;
    typedef typename Input1ImageType::ConstPointer Input1ImagePointer;
    typedef typename Input1ImageType::RegionType Input1ImageRegionType;
    typedef typename Input1ImageType::PixelType Input1ImagePixelType;

    //GM image
    typedef TInputImage2 Input2ImageType;
    typedef typename Input2ImageType::ConstPointer Input2ImagePointer;
    typedef typename Input2ImageType::RegionType Input2ImageRegionType;
    typedef typename Input2ImageType::PixelType Input2ImagePixelType;

    //WM image
    typedef TInputImage3 Input3ImageType;
    typedef typename Input3ImageType::ConstPointer Input3ImagePointer;
    typedef typename Input3ImageType::RegionType Input3ImageRegionType;
    typedef typename Input3ImageType::PixelType Input3ImagePixelType;

    typedef TOutputImage OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename OutputImageType::RegionType OutputImageRegionType;
    typedef typename OutputImageType::PixelType OutputImagePixelType;

    //PET
    void SetInput1(const TInputImage1 * imagePET);

    //GM
    void SetInput2(const TInputImage2 * imageGM);

    //WM
    void SetInput3(const TInputImage3 * imageWM);

    void SetPSF(VectorType dSigma2) {
        this->m_dVariance = dSigma2;
    };

    void SetVerbose(bool bVerbose) {
        this->m_bVerbose = bVerbose;
    };

    void SetWM(double fWM) {
        this->m_fGivenWM = fWM;
    };

    /** ImageDimension constants */
    itkStaticConstMacro(
        InputImage1Dimension, unsigned int, TInputImage1::ImageDimension);
    itkStaticConstMacro(
        InputImage2Dimension, unsigned int, TInputImage2::ImageDimension);
    itkStaticConstMacro(
        InputImage3Dimension, unsigned int, TInputImage3::ImageDimension);
    itkStaticConstMacro(
        OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro(SameDimensionCheck1,
                    (Concept::SameDimension<itkGetStaticConstMacro(InputImage1Dimension),
                     itkGetStaticConstMacro(InputImage2Dimension)>));
    itkConceptMacro(SameDimensionCheck2,
                    (Concept::SameDimension<itkGetStaticConstMacro(InputImage1Dimension),
                     itkGetStaticConstMacro(InputImage3Dimension)>));
    itkConceptMacro(SameDimensionCheck3,
                    (Concept::SameDimension<itkGetStaticConstMacro(InputImage1Dimension),
                     itkGetStaticConstMacro(OutputImageDimension)>));
    /** End concept checking */
#endif

protected:
    MullerGartnerImageFilter();

    virtual ~MullerGartnerImageFilter() {
    }

    //void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId );
    void GenerateData();

protected:

    const static int nDimension = 3;

    /*! All internal operations are performed using the type defined as InternalPixelType.*/
    typedef Input1ImagePixelType InternalPixelType;

    typedef Input1ImageType InternalImageType;


    /*! Defines step for use with MultiplyByConstantImageFilter. */
    typedef InternalPixelType StepType;

    /*! Defines image iterator. */
    typedef itk::ImageRegionConstIterator<InternalImageType> ItType;

    /*! Defines Gaussian filter type.
     * Blurs an image with a Gaussian. The variance (sigma squared) specifies the width. */
    typedef itk::DiscreteGaussianImageFilter< InternalImageType, InternalImageType > GaussianFilterType;

    /*! Defines subtract image filter type.
     * Subtracts image 2 from image 1.
     * */
    typedef itk::SubtractImageFilter<InternalImageType, InternalImageType, InternalImageType> SubtractionImageFilterType;

    /*! Defines multiply image filter type.
     * Multiplies image 1 by image 2.
     * */
    typedef itk::MultiplyImageFilter<InternalImageType, InternalImageType, InternalImageType> MultiplyImageFilterType;

    /*! Defines divide image filter type.
     * Divides image 1 by image 2.
     * */
    typedef itk::DivideImageFilter<InternalImageType, InternalImageType, InternalImageType> DivideImageFilterType;

    /*! Defines duplicate image filter type.
     * Makes a complete copy of the input. Useful when images need to be stored at intermediate
     * stages in the pipeline.
     * */
    typedef itk::ImageDuplicator<InternalImageType> ImageDuplicatorType;

    /*! Defines threshold image filter type for binary images.
     * Currently just assumes input will be binary (0 and 1) image. Could be used for fuzzy
     * images instead.
     * */
    typedef itk::BinaryThresholdImageFilter< InternalImageType, InternalImageType >
    BinaryThresholdFilterType;

    /*! Defines threshold image filter type for  images.
     * */

    typedef itk::ThresholdImageFilter< InternalImageType > ThresholdFilterType;

    /*! Defines ball structuring element for erosion.*/
    typedef itk::BinaryBallStructuringElement< InternalPixelType, nDimension > StructuringElementType;
    /*! Defines erosion filter, using ball structuring element */
    typedef itk::BinaryErodeImageFilter< InternalImageType, InternalImageType, StructuringElementType > ErodeFilterType;


private:
    MullerGartnerImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    typename GaussianFilterType::Pointer m_filterGaussian;
    typename GaussianFilterType::Pointer m_filterGaussian2;

    typename SubtractionImageFilterType::Pointer m_subFilter;

    typename MultiplyImageFilterType::Pointer m_multiplyFilter;
    typename DivideImageFilterType::Pointer m_divideFilter;

    typename MultiplyImageFilterType::Pointer m_multFilter;
    typename ImageDuplicatorType::Pointer m_duplicator;
    typename BinaryThresholdFilterType::Pointer m_thresholdFilter;
    typename ThresholdFilterType::Pointer m_thresholdFilter2;
    typename ErodeFilterType::Pointer m_binaryErodeFilter;

    typename InternalImageType::Pointer m_imgErodedWM;

    VectorType m_dVariance;

    bool m_bVerbose;
    double m_fGivenWM;

};

} // end namespace petpvc

#ifndef ITK_MANUAL_INSTANTIATION
#include "petpvcMullerGartnerImageFilter.txx"
#endif

#endif
