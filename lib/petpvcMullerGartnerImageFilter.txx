/*
   petpvcMullerGartnerImageFilter.txx

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
#ifndef __PETPVCMULLERGARTNERIMAGEFILTER_TXX
#define __PETPVCMULLERGARTNERIMAGEFILTER_TXX

#include "petpvcMullerGartnerImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace petpvc
{

/**
 * Constructor
 */
template <class TInputImage1, class TInputImage2, class TInputImage3,
         class TOutputImage>
MullerGartnerImageFilter<TInputImage1, TInputImage2, TInputImage3, TOutputImage>
::MullerGartnerImageFilter()
{
    this->SetNumberOfRequiredInputs(3);
    this->InPlaceOff();

    this->m_dVariance = 0.0;
    this->m_bVerbose = false;

    this->m_fGivenWM = 0.0;

    m_filterGaussian = GaussianFilterType::New();
    m_filterGaussian2 = GaussianFilterType::New();

    m_subFilter = SubtractionImageFilterType::New();

    m_multiplyFilter = MultiplyImageFilterType::New();
    m_divideFilter = DivideImageFilterType::New();
    m_multFilter = MultiplyImageFilterType::New();
    m_duplicator = ImageDuplicatorType::New();
    m_thresholdFilter = BinaryThresholdFilterType::New();
    m_binaryErodeFilter = ErodeFilterType::New();
    m_thresholdFilter2 = ThresholdFilterType::New();

    m_imgErodedWM = InternalImageType::New();
}

template <class TInputImage1, class TInputImage2, class TInputImage3,
         class TOutputImage>
void
MullerGartnerImageFilter<TInputImage1, TInputImage2, TInputImage3, TOutputImage>
::SetInput1(const TInputImage1 * imagePET)
{
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput(0, const_cast<TInputImage1 *> (imagePET));
}

template <class TInputImage1, class TInputImage2, class TInputImage3,
         class TOutputImage>
void
MullerGartnerImageFilter<TInputImage1, TInputImage2, TInputImage3, TOutputImage>
::SetInput2(const TInputImage2 * imageGM)
{
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput(1, const_cast<TInputImage2 *> (imageGM));
}

template <class TInputImage1, class TInputImage2, class TInputImage3,
         class TOutputImage>
void
MullerGartnerImageFilter<TInputImage1, TInputImage2, TInputImage3, TOutputImage>
::SetInput3(const TInputImage3 * imageWM)
{
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput(2, const_cast<TInputImage3 *> (imageWM));
}

/**
 * GenerateData Performs the Muller Gartner PV correction.
 */
template <class TInputImage1, class TInputImage2, class TInputImage3, class TOutputImage>
void
MullerGartnerImageFilter<TInputImage1, TInputImage2, TInputImage3, TOutputImage>
::GenerateData()
{
    //Get input images
    Input1ImagePointer inputPET
        = dynamic_cast<const TInputImage1*> (ProcessObject::GetInput(0));
    Input2ImagePointer inputGM
        = dynamic_cast<const TInputImage2*> (ProcessObject::GetInput(1));
    Input3ImagePointer inputWM
        = dynamic_cast<const TInputImage3*> (ProcessObject::GetInput(2));
    OutputImagePointer outputPtr = this->GetOutput(0);

    //Setup blurring filters
    m_filterGaussian->SetVariance(this->m_dVariance);
    m_filterGaussian2->SetVariance(this->m_dVariance);


	//Pick size of erosion depending on voxel size.
	typename Input1ImageType::SpacingType voxelSize = inputPET->GetSpacing();
	float avgDim = (voxelSize[0]+ voxelSize[1] + voxelSize[2] )/3;
	float elementRadius;

	if ( avgDim <= 1.0 )
		elementRadius = 6;
	else
		if ( avgDim >= 2.5 )
			elementRadius = 3;
		else
			elementRadius = floor(( avgDim * -2.0 + 8.0 )+0.5f);

    if (this->m_bVerbose) {
        std::cout << "Structuring element : " << elementRadius << std::endl;
    }

    //Setup WM erosion filter
    StructuringElementType structuringElement;
    structuringElement.SetRadius( elementRadius );
    structuringElement.CreateStructuringElement();
    m_binaryErodeFilter->SetKernel(structuringElement);
    m_binaryErodeFilter->SetInput(inputWM);
    m_binaryErodeFilter->SetErodeValue(1);

    //Threshold erosion filter to only include values = 1
    m_thresholdFilter->SetInput(m_binaryErodeFilter->GetOutput());
    m_thresholdFilter->SetLowerThreshold(1);
    m_thresholdFilter->SetUpperThreshold(1);
    m_thresholdFilter->SetOutsideValue(0);
    m_thresholdFilter->SetInsideValue(1);
    m_thresholdFilter->Update();

    m_imgErodedWM = m_thresholdFilter->GetOutput();

    //Multiply PET by eroded WM
    m_multiplyFilter->SetInput1(inputPET);
    m_multiplyFilter->SetInput2(m_imgErodedWM);
    m_multiplyFilter->Update();

    //Calculate sum of previous multiplication and sum of the eroded WM mask.
    ItType iteratorVirtualWM(m_multiplyFilter->GetOutput(), m_multiplyFilter->GetOutput()->GetRequestedRegion());
    ItType iteratorErodedWM(m_thresholdFilter->GetOutput(), m_thresholdFilter->GetOutput()->GetRequestedRegion());

    double sumVWM = 0.0;
    double sumErodedWM = 0.0;

    while (!iteratorVirtualWM.IsAtEnd() ) {
        sumVWM += iteratorVirtualWM.Value();
        ++iteratorVirtualWM;
    }

    while (!iteratorErodedWM.IsAtEnd() ) {
        sumErodedWM += iteratorErodedWM.Value();
        ++iteratorErodedWM;
    }

    //Calculate mean value in WM.
    double fWhiteMatterMeanValue = sumVWM / sumErodedWM;

    if (this->m_fGivenWM != 0.0) {
        fWhiteMatterMeanValue = this->m_fGivenWM;
    }

    if (this->m_bVerbose) {
        std::cout << "White Matter mean value : " << fWhiteMatterMeanValue << std::endl;
    }

    //Apply WM value to WM mask.
    m_multFilter->SetConstant(fWhiteMatterMeanValue);
    m_multFilter->SetInput(inputWM);

    //Smooth output of previous step by resolution of PET.
    m_filterGaussian->SetInput(m_multFilter->GetOutput());

    //Remove Smooth WM estimate from PET
    m_subFilter->SetInput1(inputPET);
    m_subFilter->SetInput2(m_filterGaussian->GetOutput());

    //Clip to GM
    m_multiplyFilter->SetInput1(m_subFilter->GetOutput());
    m_multiplyFilter->SetInput2(inputGM);

    //Smooth GM mask
    m_filterGaussian2->SetInput(inputGM);

    //Divide vGM by smoothed GM mask (apply GM correction factors)
    m_divideFilter->SetInput1(m_multiplyFilter->GetOutput());
    m_divideFilter->SetInput2(m_filterGaussian2->GetOutput());
    m_divideFilter->Update();

    //Clean up areas where you get infs due to divide by zeros.
    m_thresholdFilter2->SetInput(m_divideFilter->GetOutput());
    m_thresholdFilter2->ThresholdOutside(0, fWhiteMatterMeanValue * 10.0);
    m_thresholdFilter2->SetOutsideValue(0);
    m_thresholdFilter2->Update();

    this->GraftOutput(m_thresholdFilter2->GetOutput());

}

} // end namespace petpvc

#endif
