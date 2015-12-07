/*
   petpvcVanCittertPVCImageFilter.txx

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

#ifndef __PETPVCVANCITTERTPVCIMAGEFILTER_TXX
#define __PETPVCVANCITTERTPVCIMAGEFILTER_TXX

#include "petpvcVanCittertPVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage >
VanCittertPVCImageFilter< TInputImage >
::VanCittertPVCImageFilter()
{
    this->m_nIterations = 30;
    this->m_bVerbose = false;
    this->m_fAlpha = 1.5;
    this->m_fStopCriterion = 0.01;
}

template< class TInputImage >
void VanCittertPVCImageFilter< TInputImage >
::GenerateData()
{
    typename TInputImage::ConstPointer input = this->GetInput();
    typename TInputImage::Pointer output = this->GetOutput();

    InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));

    //Gaussian smoothing
    typename BlurringFilterType::Pointer blurFilter = BlurringFilterType::New();
    typename BlurringFilterType::Pointer blurFilter2 = BlurringFilterType::New();
    //Multiplies two images together.
    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
    typename SubFilterType::Pointer subFilter = SubFilterType::New();
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();

    duplicator->SetInputImage( pPET );
    duplicator->Update();

    typename TInputImage::Pointer imageEstimate;
    //Set image estimate to the original PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    typename TInputImage::Pointer imagePrev;
    //Set imagePrev to the original PET data for the first iteration.
    imagePrev = duplicator->GetOutput();

    blurFilter->SetVariance( this->GetPSF() );
    blurFilter2->SetVariance( this->GetPSF() );

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
            std::cout << n << "\t" << fCurrentEval << std::endl;
            n++;

            if ( fCurrentEval < this->m_fStopCriterion )
                bStopped = true;
    }
    std::cout << std::endl;


  
    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageEstimate.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
