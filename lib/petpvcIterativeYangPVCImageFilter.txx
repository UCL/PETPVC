/*
   petpvcIterativeYangPVCImageFilter.txx

   Authors:     Benjamin A. Thomas
                Kris Thielemans (minor modifications)

   Copyright 2013-2014 Institute of Nuclear Medicine, University College London.

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

#ifndef __PETPVCITERATIVEYANGPVCIMAGEFILTER_TXX
#define __PETPVCITERATIVEYANGPVCIMAGEFILTER_TXX

#include "petpvcIterativeYangPVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
IterativeYangPVCImageFilter< TInputImage, TMaskImage>
::IterativeYangPVCImageFilter()
{
    this->m_nIterations = 10;
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void IterativeYangPVCImageFilter< TInputImage, TMaskImage>
::GenerateData()
{

    itk::ImageToImageFilterCommon::SetGlobalDefaultCoordinateTolerance( 1e-2 );
    itk::ImageToImageFilterCommon::SetGlobalDefaultDirectionTolerance( 1e-2 );
    
    typename TInputImage::ConstPointer input = this->GetInput();
    typename TInputImage::Pointer output = this->GetOutput();

    typename FuzzyCorrFilterType::Pointer pFuzzyCorrFilter = FuzzyCorrFilterType::New();

    InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));
    MaskImagePointer pMask = dynamic_cast<const TMaskImage*> (ProcessObject::GetInput(1));

    pFuzzyCorrFilter->SetInput( pMask );

    //Calculate Fuzziness.
    if ( this->m_bVerbose ) {
      std::cout << "Start fuzziness calculation" << std::endl;
    }

    pFuzzyCorrFilter->Update();

    if ( this->m_bVerbose ) {
      std::cout << "matrix:\n" << pFuzzyCorrFilter->GetMatrix() << std::endl;
    }

    //Get fuzziness correction factors.
    vnl_matrix<float> matFuzzyCorr = pFuzzyCorrFilter->GetMatrix();
    vnl_vector<float> vecRegSize = pFuzzyCorrFilter->GetSumOfRegions();

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

    typename TInputImage::Pointer imageExtractedRegion;// = InputImagePointer::New();

    float fSumOfPETReg;

    //Vector to contain the current estimate of the regional mean values.
    vnl_vector<float> vecRegMeansCurrent;
    vecRegMeansCurrent.set_size(nClasses);

    //Vector to contain the estimated means after GTM correction.
    vnl_vector<float> vecRegMeansUpdated;
    vecRegMeansUpdated.set_size(nClasses);


    //Applying the Yang correction step:

    typename TInputImage::Pointer imageYang;
    typename TInputImage::Pointer imageEstimate;
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

    typename MultiplyFilterType::Pointer multiplyFilter2 = MultiplyFilterType::New();
    typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
    typename BlurringFilterType::Pointer pBlurFilter = BlurringFilterType::New();


    duplicator->SetInputImage( pPET );
    duplicator->Update();

    //Set image estimate to the original PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    int nNumOfIters =  this->m_nIterations;

    for (int k = 1; k <= nNumOfIters; k++) {

        if ( this->m_bVerbose ) {
            if (k == 1) {
                std::cout << std::endl << "Iteration:  " << std::endl;
            }

            std::cout << k << "  "; // << std::endl;
            std::flush(std::cout);
        }

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
            imageExtractedRegion->SetDirection(imageEstimate->GetDirection());
            imageExtractedRegion->UpdateOutputData();


            //Multiply current image estimate by region mask. To clip PET values
            //to mask.
            multiplyFilter->SetInput1(imageEstimate);
            multiplyFilter->SetInput2(imageExtractedRegion);

            statsFilter->SetInput(multiplyFilter->GetOutput());
            statsFilter->Update();

            //Get sum of the clipped image.
            fSumOfPETReg = statsFilter->GetSum();

            //Place regional mean into vector.
			float fNewRegMean = std::max( (float) (fSumOfPETReg / vecRegSize.get(i - 1)), (float)0.0);
            vecRegMeansCurrent.put(i - 1, fNewRegMean );

            //std::cout << std::endl << "Sum = " << fSumOfPETReg << " , " << "Mean = " << vecRegMeansCurrent.get( i-1 ) << " , Size = " << vecRegSize.get(i - 1) << std::endl;

        }


        //Apply fuzziness correction to current mean value estimates.
        vecRegMeansUpdated = vnl_matrix_inverse<float>( matFuzzyCorr )
                             * vecRegMeansCurrent;

        //std::cout << vecRegMeansCurrent << std::endl;
        if ( this->m_bVerbose ) {
            std::cout << vecRegMeansUpdated << std::endl;
        }
        /*
        for (int n = 0; n < vecRegMeansUpdated.size(); n++) {
        	float fNewRegMean = fmax ( vecRegMeansUpdated.get(n) , 0.0 );
        	vecRegMeansUpdated.put(n, fNewRegMean);
        }*/

        //std::cout << vecRegMeansUpdated << std::endl;

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

            //Multiply current image estimate by region mask. To clip PET values
            //to mask.
            multiplyFilter->SetInput1( vecRegMeansUpdated.get(i-1) );
            multiplyFilter->SetInput2( imageExtractedRegion );
            multiplyFilter->Update();

            //If this is the first region, create imageYang,
            //else add the current region to the previous contents of imageYang.
            if (i == 1) {
                imageYang = multiplyFilter->GetOutput();
                imageYang->DisconnectPipeline();
            } else {
                addFilter->SetInput1(imageYang);
                addFilter->SetInput2(multiplyFilter->GetOutput());
                addFilter->Update();

                imageYang = addFilter->GetOutput();
            }

        }

        //Takes the original PET data and the pseudo PET image, calculates the
        //correction factors  and returns the PV-corrected PET image.

        pBlurFilter->SetInput(imageYang);
        pBlurFilter->SetVariance( this->GetPSF() );

        //Take ratio of pseudo PET and smoothed pseudo PET. These are the correction
        //factors.
        divideFilter->SetInput1( imageYang );
        divideFilter->SetInput2(pBlurFilter->GetOutput());

        //Multiply original PET by correction factors.
        multiplyFilter2->SetInput1( pPET );
        multiplyFilter2->SetInput2(divideFilter->GetOutput());
        multiplyFilter2->Update();


        imageEstimate = multiplyFilter2->GetOutput();
        imageEstimate->UpdateOutputData();
        imageEstimate->DisconnectPipeline();
    }

    if ( this->m_bVerbose ) {
        std::cout << std::endl;
    }

    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageEstimate.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
