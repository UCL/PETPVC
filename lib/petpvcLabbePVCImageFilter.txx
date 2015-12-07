/*
   petpvcLabbePVCImageFilter.txx

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

#ifndef __PETPVCLABBEPVCIMAGEFILTER_TXX
#define __PETPVCLABBEPVCIMAGEFILTER_TXX

#include "petpvcLabbePVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
LabbePVCImageFilter< TInputImage, TMaskImage>
::LabbePVCImageFilter()
{
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void LabbePVCImageFilter< TInputImage, TMaskImage>
::GenerateData()
{
    typename TInputImage::ConstPointer input = this->GetInput();
    typename TInputImage::Pointer output = this->GetOutput();

    typename LabbeImageFilterType::Pointer pLabbe = LabbeImageFilterType::New();

    InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));
    MaskImagePointer pMask = dynamic_cast<const TMaskImage*> (ProcessObject::GetInput(1));

    pLabbe->SetInput( pMask );
    pLabbe->SetPSF( this->GetPSF() );
    //Calculate Labbe.
    try {
        pLabbe->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tCannot calculate Labbe"
                  << std::endl;
    }

    if ( this->m_bVerbose ) {
        std::cout << pLabbe->GetMatrix() << std::endl;
    }

    /////////////////////////////////////////////

    //Get mask image size.
    typename MaskImageType::SizeType imageSize =
        pMask->GetLargestPossibleRegion().GetSize();

    int nClasses = 0;

    //If mask is not 4D, then quit.
    if (imageSize.Dimension == 4) {
        nClasses = imageSize[3];
    } else {
        std::cerr << "[Error]\tMask file must be 4-D!"
                  << std::endl;
    }

    MaskSizeType desiredStart;
    desiredStart.Fill(0);
    MaskSizeType desiredSize = imageSize;

    //Extract filter used to extract 3D volume from 4D file.
    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput( pMask );
    extractFilter->SetDirectionCollapseToIdentity(); // This is required.

    //Stats. filter used to calculate statistics for an image.
    typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();

    //Multiplies two images together.
    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();

	//Smooth the pseudo PET by the PSF.
    typename BlurringFilterType::Pointer blurFilter = BlurringFilterType::New();

    blurFilter->SetVariance( this->GetPSF() );

    typename TInputImage::Pointer imageExtractedRegion;// = InputImagePointer::New();

    float fSumOfPETReg;

    //Vector to contain the current estimate of the regional mean values.
    vnl_vector<float> vecRegMeansCurrent;
    vecRegMeansCurrent.set_size(nClasses);

    //Vector to contain the estimated means after Labbe correction.
    vnl_vector<float> vecRegMeansUpdated;
    vecRegMeansUpdated.set_size(nClasses);

    for (int i = 1; i <= nClasses; i++) {

        //Starts reading from 4D volume at index (0,0,0,i) through to
        //(maxX, maxY, maxZ,0), i.e. one 3D brain mask.
        desiredStart[3] = i - 1;
        desiredSize[3] = 0;

        //Get region mask.
        MaskRegionType maskReg;
        maskReg.SetSize(desiredSize );
        maskReg.SetIndex(0,desiredStart[0] );
        maskReg.SetIndex(1,desiredStart[1] );
        maskReg.SetIndex(2,desiredStart[2] );
        maskReg.SetIndex(3,desiredStart[3] );

        extractFilter->SetExtractionRegion( maskReg );
        extractFilter->Update();

        imageExtractedRegion = extractFilter->GetOutput();
        imageExtractedRegion->SetDirection( pPET->GetDirection() );
        imageExtractedRegion->UpdateOutputData();

		blurFilter->SetInput( imageExtractedRegion );

        //Multiply current image estimate by smoothed region mask. To clip PET values
        //to mask.
        multiplyFilter->SetInput1( pPET );
        multiplyFilter->SetInput2( blurFilter->GetOutput() );

        statsFilter->SetInput(multiplyFilter->GetOutput());
        statsFilter->Update();

        //Get sum of the clipped image.
        fSumOfPETReg = statsFilter->GetSum();

        //Place regional mean into vector.
        vecRegMeansCurrent.put(i - 1, fSumOfPETReg / pLabbe->GetSumOfRegions().get(i - 1));
        //std::cout << "Sum = " << fSumOfPETReg << " , Mean = " << vecRegMeansCurrent.get(i-1) << " Total vox. = " << gtmFilter->GetSumOfRegions().get( i-1 ) << std::endl;

    }

    //Apply Labbe to regional mean values.
    vecRegMeansUpdated = vnl_matrix_inverse<float>(pLabbe->GetMatrix()) * vecRegMeansCurrent;

    if ( this->m_bVerbose ) {
        std::cout << std::endl << "Regional means:" << std::endl;
        std::cout << vecRegMeansCurrent << std::endl << std::endl;

        std::cout << "Labbe:" << std::endl;
        pLabbe->GetMatrix().print(std::cout);

        std::cout << std::endl << "Corrected means:" << std::endl;

    }

    std::cout << vecRegMeansUpdated << std::endl;
	this->m_vecRegMeansPVCorr = vecRegMeansUpdated;

    this->AllocateOutputs();

    ImageAlgorithm::Copy(input.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                         output->GetRequestedRegion() );


}

}// end namespace


#endif
