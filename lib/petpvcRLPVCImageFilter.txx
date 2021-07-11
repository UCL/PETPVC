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
float RichardsonLucyPVCImageFilter< TInputImage >
::GetZeroThreshold( typename TInputImage::ConstPointer img )
{
    // Default threshold to zero values at:
    const float fZeroThreshold = 1e-4f;

    // Calculate image statistics
    typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
    statsFilter->SetInput( img );
    statsFilter->Update();
    
    // Return default or (image max * default ) as new threshold
    const float fNewThreshold = std::min( fZeroThreshold, statsFilter->GetMaximum() * fZeroThreshold );

    // Set to zero if somehow new threshold is negative!
    return std::max( fNewThreshold , 0.0f );
}

template< class TInputImage >
void RichardsonLucyPVCImageFilter< TInputImage >
::GenerateData()
{
    this->SetGlobalDefaultCoordinateTolerance( 1e-2 );
    this->SetGlobalDefaultDirectionTolerance( 1e-2 );

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
    imageEstimate->DisconnectPipeline();

    // Create a blank image with same dimensions as input to use for division
    duplicator->SetInputImage(imageEstimate);
    duplicator->Update();

    typedef itk::CastImageFilter< TInputImage, TInputImage > CastImageFilterType;
    typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
    castFilter->SetInput( duplicator->GetOutput() );
    castFilter->Update();
    typename TInputImage::Pointer blankImage = castFilter->GetOutput();

    typedef typename TInputImage::PixelType PixelType;
    blankImage->FillBuffer( itk::NumericTraits< PixelType >::Zero );

    // Image to store result of division
    typename TInputImage::Pointer dividedImage;

    int nMaxNumOfIters =  this->m_nIterations;
    int n=1;

    bool bStopped = false;
	
    while ( ( n <= nMaxNumOfIters ) && ( !bStopped ) ) {
         float fLog = 0.0;  
            // f_k * h
            blurFilter->SetInput( imageEstimate );
            blurFilter->Update();

            // Perform f(x) / [ f_k(x) * h ] voxel-by-voxel as itk::DivideFilterType sets image
            // to max if denominator is 0.
            ConstImageIterator numeratorIt( thresholdFilter->GetOutput(),  thresholdFilter->GetOutput()->GetLargestPossibleRegion() );
            ConstImageIterator denomIt( blurFilter->GetOutput(),  blurFilter->GetOutput()->GetLargestPossibleRegion() );

            duplicator->SetInputImage(blankImage);
            duplicator->Update();
            dividedImage = duplicator->GetOutput();
            dividedImage->DisconnectPipeline();

            ImageRegionIterator<TInputImage> divIt( dividedImage, dividedImage->GetLargestPossibleRegion() );

            numeratorIt.GoToBegin();
            denomIt.GoToBegin();
            divIt.GoToBegin();

            // Get zero threshold
            const float fSmallNum = this->GetZeroThreshold( blurFilter->GetOutput() );

            while ( !divIt.IsAtEnd() )
    	    {
				if ( denomIt.Get() > fSmallNum )
                    divIt.Set(numeratorIt.Get() / denomIt.Get());
                else
                    divIt.Set(itk::NumericTraits< PixelType >::Zero);
			
        		++numeratorIt;
                ++denomIt;
                ++divIt;
        	}

            // Reblur correction factors            
            blurFilter2->SetInput( dividedImage );
            blurFilter2->Update();

            // Multiply current image estimate by reblurred correction factors
			multiplyFilter->SetInput1( imageEstimate );
			multiplyFilter->SetInput2( blurFilter2->GetOutput() );
            multiplyFilter->Update();
            
            // Update image estimate 
            imageEstimate = multiplyFilter->GetOutput();
            imageEstimate->DisconnectPipeline();

            ConstImageIterator origIt( thresholdFilter->GetOutput(),  thresholdFilter->GetOutput()->GetLargestPossibleRegion() );
            ConstImageIterator currIt( imageEstimate, imageEstimate->GetLargestPossibleRegion() );

            origIt.GoToBegin();
            currIt.GoToBegin();

            while ( !currIt.IsAtEnd() )
    	    {
				if (currIt.Get() > 0.0) {
                    fLog += origIt.Get() * log(currIt.Get()) - currIt.Get();
				}

        		++currIt;
                ++origIt;
        	}

            float fCurrentEval = fLog;
            std::cout << n << "\t" << fCurrentEval << std::endl;         
			n++;
    }
    std::cout << std::endl;

    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageEstimate.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
