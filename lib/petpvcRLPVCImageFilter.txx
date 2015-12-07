/*
   petpvcRLPVCImageFilter.txx

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

#ifndef __PETPVCRLPVCIMAGEFILTER_TXX
#define __PETPVCRLPVCIMAGEFILTER_TXX

#include "petpvcRLPVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage >
RichardsonLucyPVCImageFilter< TInputImage >
::RichardsonLucyPVCImageFilter()
{
    this->m_nIterations = 10;
    this->m_bVerbose = false;
    this->m_fStopCriterion = -3e+06;
}

template< class TInputImage >
void RichardsonLucyPVCImageFilter< TInputImage >
::GenerateData()
{
    typename TInputImage::ConstPointer input = this->GetInput();
    typename TInputImage::Pointer output = this->GetOutput();

    InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));

    //Gaussian smoothing
    typename BlurringFilterType::Pointer blurFilter = BlurringFilterType::New();
	typename BlurringFilterType::Pointer blurFilter2 = BlurringFilterType::New();

    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
	typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();

    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
    typename SubFilterType::Pointer subFilter = SubFilterType::New();
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();

    blurFilter->SetVariance( this->GetPSF() );
	blurFilter2->SetVariance( this->GetPSF() );

    thresholdFilter->ThresholdBelow( 0 );
    thresholdFilter->SetOutsideValue( 0 );
	thresholdFilter->SetInput( pPET );  
	thresholdFilter->Update();
	duplicator->SetInputImage( thresholdFilter->GetOutput() );
    duplicator->Update();

    typename TInputImage::Pointer imageEstimate;
    //Set image estimate to the original non-negative PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    int nMaxNumOfIters =  this->m_nIterations;
    int n=1;

    bool bStopped = false;

	
    while ( ( n <= nMaxNumOfIters ) && ( !bStopped ) ) {
         float fLog = 0.0;   
            blurFilter->SetInput( imageEstimate );
			divideFilter->SetInput1( thresholdFilter->GetOutput() );
			divideFilter->SetInput2( blurFilter->GetOutput() );
			multiplyFilter->SetInput1( imageEstimate );
			multiplyFilter->SetInput2( divideFilter->GetOutput() );
            multiplyFilter->Update();

            imageEstimate = multiplyFilter->GetOutput();
            imageEstimate->DisconnectPipeline();
            
			blurFilter2->SetInput( imageEstimate );
			blurFilter2->Update();

            ConstImageIterator origIt( thresholdFilter->GetOutput(),  thresholdFilter->GetOutput()->GetLargestPossibleRegion() );
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

            float fCurrentEval = fLog;
            std::cout << n << "\t" << fCurrentEval << std::endl;         
			n++;

            //if ( fCurrentEval < this->m_fStopCriterion )
            //   bStopped = true;
    }
    std::cout << std::endl;

    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageEstimate.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
