/*
   petpvcIntraRegVCImageFilter.txx

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

#ifndef __PETPVCINTRAREGVCIMAGEFILTER_TXX
#define __PETPVCINTRAREGVCIMAGEFILTER_TXX

#include "petpvcIntraRegVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
IntraRegVCImageFilter< TInputImage, TMaskImage>
::IntraRegVCImageFilter()
{
    this->m_nIterations = 10;
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void IntraRegVCImageFilter< TInputImage, TMaskImage>
::GenerateData()
{
    typename TInputImage::ConstPointer input = this->GetInput();
    typename TInputImage::Pointer output = this->GetOutput();

    InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));
    MaskImagePointer pMask = dynamic_cast<const TMaskImage*> (ProcessObject::GetInput(1));

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

    typename TInputImage::Pointer imageExtractedRegion;// = InputImagePointer::New();

    //Applying the Intra-regional reblurred Van-Cittert:
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
    typename SubFilterType::Pointer subFilter = SubFilterType::New();
    typename MultiplyFilterType::Pointer multiplyFilter2 = MultiplyFilterType::New();
    typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
    typename IntraRegBlurFilterType::Pointer blurFilter = IntraRegBlurFilterType::New();
	typename IntraRegBlurFilterType::Pointer blurFilter2 = IntraRegBlurFilterType::New();
    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();

    duplicator->SetInputImage( pPET );
    duplicator->Update();

    typename TInputImage::Pointer imageEstimate;
    //Set image estimate to the original PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    typename TInputImage::Pointer imagePrev;
    //Set imagePrev to the original PET data for the first iteration.
    imagePrev = duplicator->GetOutput();

    blurFilter->SetPSF( this->GetPSF() );
    blurFilter2->SetPSF( this->GetPSF() );

    thresholdFilter->ThresholdBelow( 0 );
    thresholdFilter->SetOutsideValue( 0 );

    float fSumOfPETsq = 0.0;

	ConstImageIterator inputIt( pPET, pPET->GetLargestPossibleRegion() );
    inputIt.GoToBegin();

	while ( !inputIt.IsAtEnd() )
	{
        fSumOfPETsq += inputIt.Get()*inputIt.Get();
		++inputIt;
	}

    int nMaxNumOfIters =  this->m_nIterations;
    int n=1;

    bool bStopped = false;

    int nNumOfIters =  this->m_nIterations;

	for (int i = 1; i <= nClasses; i++) {

		std::cout << "Region " << i << " : ";


		//Starts reading from 4D volume at index (0,0,0,i) through to
		//(maxX, maxY, maxZ,0), i.e. one 3D brain mask.
		desiredStart[3] = i - 1;
		desiredSize[3] = 0;

		//Get region mask.
		MaskRegionType maskReg;
		maskReg.SetSize( desiredSize );
		maskReg.SetIndex(0,desiredStart[0] );
		maskReg.SetIndex(1,desiredStart[1] );
		maskReg.SetIndex(2,desiredStart[2] );
		maskReg.SetIndex(3,desiredStart[3] );

		extractFilter->SetExtractionRegion( maskReg );
		extractFilter->Update();

		imageExtractedRegion = extractFilter->GetOutput();
		imageExtractedRegion->SetDirection(imageEstimate->GetDirection());
		imageExtractedRegion->UpdateOutputData();

		blurFilter->SetMaskInput( imageExtractedRegion );
		blurFilter2->SetMaskInput( imageExtractedRegion );
	

    	while ( ( n <= nMaxNumOfIters ) && ( !bStopped ) ) {
            
            float fSumOfDiffsq = 0.0;

            duplicator->SetInputImage( imageEstimate );
            duplicator->Update();
            imagePrev = duplicator->GetOutput();

            blurFilter->SetInput( imageEstimate );
            subFilter->SetInput1( pPET );
            subFilter->SetInput2( blurFilter->GetOutput() );
            blurFilter2->SetInput( subFilter->GetOutput() ); 
            multiplyFilter->SetConstant( this->m_fAlpha );
            multiplyFilter->SetInput( blurFilter2->GetOutput() );
            addFilter->SetInput1( imageEstimate );
            addFilter->SetInput2( multiplyFilter->GetOutput() );
            thresholdFilter->SetInput( addFilter->GetOutput() );  
            thresholdFilter->Update();

            imageEstimate = thresholdFilter->GetOutput();
            imageEstimate->DisconnectPipeline();
            
            ConstImageIterator currIt( imageEstimate, imageEstimate->GetLargestPossibleRegion() );
            ConstImageIterator prevIt( imagePrev, imagePrev->GetLargestPossibleRegion() );

            currIt.GoToBegin();
            prevIt.GoToBegin();

            while ( !currIt.IsAtEnd() )
    	    {
                float diff = currIt.Get() - prevIt.Get();
                fSumOfDiffsq += diff*diff;
        		++currIt;
                ++prevIt;
        	}

            float fCurrentEval = sqrt( fSumOfDiffsq ) / sqrt( fSumOfPETsq );
            //std::cout << n << "\t" << fCurrentEval << std::endl;
		    std::cout << n << " " << std::flush; 
            n++;

            if ( fCurrentEval < this->m_fStopCriterion )
                bStopped = true;
    	}
		
		n=1;
		bStopped = false;
        std::cout << std::endl;
		       
		

	}// End of regions


    if ( this->m_bVerbose ) {
        std::cout << std::endl;
    }

    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageEstimate.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
