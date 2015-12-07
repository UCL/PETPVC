/*
   petpvcIntraRegRLImageFilter.txx

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

#ifndef __PETPVCINTRAREGRLIMAGEFILTER_TXX
#define __PETPVCINTRAREGRLIMAGEFILTER_TXX

#include "petpvcIntraRegRLImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
IntraRegRLImageFilter< TInputImage, TMaskImage>
::IntraRegRLImageFilter()
{
    this->m_nIterations = 10;
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void IntraRegRLImageFilter< TInputImage, TMaskImage>
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

    //Applying the Intra-regional Richardson-Lucy:
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
	typename AddFilterType::Pointer addFilter2 = AddFilterType::New();
    typename SubFilterType::Pointer subFilter = SubFilterType::New();
    typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
    typename IntraRegBlurFilterType::Pointer blurFilter = IntraRegBlurFilterType::New();
	typename IntraRegBlurFilterType::Pointer blurFilter2 = IntraRegBlurFilterType::New();
    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    typename MultiplyFilterType::Pointer multiplyFilter2 = MultiplyFilterType::New();
	typename MultiplyFilterType::Pointer clipToRegionFilter = MultiplyFilterType::New();
	typename MultiplyFilterType::Pointer clipToRegionFilter2 = MultiplyFilterType::New();
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();

    thresholdFilter->ThresholdBelow( 0 );
    thresholdFilter->SetOutsideValue( 0 );
	thresholdFilter->SetInput( pPET );  
	thresholdFilter->Update();
	duplicator->SetInputImage( thresholdFilter->GetOutput() );
    duplicator->Update();

    typename TInputImage::Pointer imageThresholded;
    //Set image to the original PET data after non-negativity constraint.
    imageThresholded = duplicator->GetOutput();

    typename TInputImage::Pointer imageEstimate;
    //Set image estimate to the original PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    typename TInputImage::Pointer imageOutput;
    //Create image to hold output.
    imageOutput = duplicator->GetOutput();

    blurFilter->SetPSF( this->GetPSF() );
    blurFilter2->SetPSF( this->GetPSF() );

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

		clipToRegionFilter->SetInput1( imageThresholded );
		clipToRegionFilter->SetInput2( imageExtractedRegion );

		clipToRegionFilter2->SetInput1( imageThresholded );
		clipToRegionFilter2->SetInput2( imageExtractedRegion );
		
		imageEstimate= clipToRegionFilter2->GetOutput();
	
    while ( ( n <= nMaxNumOfIters ) ) {
            float fLog = 0.0;  
            blurFilter->SetInput( imageEstimate );
			divideFilter->SetInput1( clipToRegionFilter->GetOutput() );
			divideFilter->SetInput2( blurFilter->GetOutput() );
			multiplyFilter->SetInput1( imageEstimate );
			multiplyFilter->SetInput2( divideFilter->GetOutput() );
            multiplyFilter->Update();

            imageEstimate = multiplyFilter->GetOutput();
            imageEstimate->DisconnectPipeline();
            
			blurFilter2->SetInput( imageEstimate );
			blurFilter2->Update();

            ConstImageIterator origIt( clipToRegionFilter->GetOutput(),  clipToRegionFilter->GetOutput()->GetLargestPossibleRegion() );
            ConstImageIterator currIt(  blurFilter2->GetOutput(), blurFilter2->GetOutput()->GetLargestPossibleRegion() );

            origIt.GoToBegin();
            currIt.GoToBegin();

            while ( !currIt.IsAtEnd() )
    	    {
				if (currIt.Get() > 0.0) {
                	float diff = origIt.Get() * log(currIt.Get()) - currIt.Get();
	                fLog += diff;
				}
        		++currIt;
                ++origIt;
        	}

            //float fCurrentEval = fLog;
            //std::cout << n << "\t" << fCurrentEval << std::endl;
			std::cout << n << " " << std::flush;         
			n++;

    }

			//If this is the first region, create imageOutput,
            //else add the current region to the previous contents of imageOutput.
            if (i == 1) {
                imageOutput = imageEstimate;
                imageOutput->DisconnectPipeline();
            } else {
                addFilter2->SetInput1( imageOutput );
                addFilter2->SetInput2( imageEstimate );
                addFilter2->Update();

                imageOutput = addFilter2->GetOutput();
            }
		
		n=1;

		std::cout << std::endl;

	}// End of regions


    if ( this->m_bVerbose ) {
        std::cout << std::endl;
    }

    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageOutput.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
