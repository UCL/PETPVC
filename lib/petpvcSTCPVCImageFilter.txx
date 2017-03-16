/*
   petpvcSTCPVCImageFilter.txx

   Author:      Benjamin A. Thomas

   Copyright 2017 Institute of Nuclear Medicine, University College London.

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

#ifndef __PETPVCSTCPVCIMAGEFILTER_TXX
#define __PETPVCSTCPVCIMAGEFILTER_TXX

#include "petpvcSTCPVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
STCPVCImageFilter< TInputImage, TMaskImage>
::STCPVCImageFilter()
{
    this->m_nIterations = 10;
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void STCPVCImageFilter< TInputImage, TMaskImage>
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

    //If mask is not 3D, then quit.
    if (imageSize.Dimension != 3) {
        std::cerr << "[Error]\tMask file must be 3-D!"
                  << std::endl;
    }

    MaskSizeType desiredStart;
    desiredStart.Fill(0);
    MaskSizeType desiredSize = imageSize;

    //Stats. filter used to calculate statistics for an image.
    typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();

    //Multiplies two images together.
    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();

    typename TInputImage::Pointer imageExtractedRegion;

    typename LabelStatisticsFilterType::Pointer labelStatsFilter = LabelStatisticsFilterType::New();
    typename BinaryThresholdImageFilterType::Pointer binThreshFilter = BinaryThresholdImageFilterType::New();

    binThreshFilter->SetInput( pMask );
    binThreshFilter->SetInsideValue(1);
    binThreshFilter->SetOutsideValue(0);

    labelStatsFilter->SetLabelInput( pMask );
    labelStatsFilter->SetInput( pPET );
    labelStatsFilter->Update();

    int numOfLabels = labelStatsFilter->GetNumberOfLabels();
    nClasses = numOfLabels;

    if ( this->m_bVerbose )
        std::cout << "Number of labels: " << nClasses << std::endl;

    float fSumOfPETReg;

    //Vector to contain the current estimate of the regional mean values.

    vnl_vector<float> vecRegMeansCurrent(nClasses);
    vecRegMeansCurrent.fill(0);

    //Vector to contain the new estimated means.
    vnl_vector<float> vecRegMeansUpdated(nClasses);
    vecRegMeansUpdated.fill(0);

    //Vector to hold size of regions.
    vnl_vector<float> vecRegSize(nClasses);
    vecRegSize.fill(0);

    //For applying the STC correction:

    typename TInputImage::Pointer imageRec;
    typename TInputImage::Pointer imageBackground;

    typename TInputImage::Pointer imageEstimate;
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
    typename SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

    typename MultiplyFilterType::Pointer multiplyFilter2 = MultiplyFilterType::New();
    typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
    typename BlurringFilterType::Pointer pBlurFilter = BlurringFilterType::New();

    duplicator->SetInputImage( pPET );
    duplicator->Update();

    //Set image estimate to the original PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    pBlurFilter->SetVariance( this->GetPSF() );

    int i=0;
    //Calculate recovery factors

    for( typename ValidLabelValuesType::const_iterator vIt=labelStatsFilter->GetValidLabelValues().begin();
            vIt != labelStatsFilter->GetValidLabelValues().end(); ++vIt)
    {
        if ( labelStatsFilter->HasLabel(*vIt) )
        {
            LabelPixelType labelValue = *vIt;
            binThreshFilter->SetLowerThreshold( labelValue );
            binThreshFilter->SetUpperThreshold( labelValue );
            binThreshFilter->Update();

            imageExtractedRegion = binThreshFilter->GetOutput();
            imageExtractedRegion->SetDirection( pPET->GetDirection() );
            imageExtractedRegion->UpdateOutputData();

            pBlurFilter->SetInput(imageExtractedRegion);
            pBlurFilter->Update();

            multiplyFilter->SetInput1( imageExtractedRegion );
            multiplyFilter->SetInput2( pBlurFilter->GetOutput() );
            multiplyFilter->Update();

            //If this is the first region, create recovery factors image,
            //else add the factors to the previous contents of imageRec.
            if (i == 0) {
                imageRec = multiplyFilter->GetOutput();
            } else {
                addFilter->SetInput1(imageRec);
                addFilter->SetInput2(multiplyFilter->GetOutput());
                addFilter->Update();

                imageRec = addFilter->GetOutput();
            }
            imageRec->DisconnectPipeline();

            i++;
        }
    }

    int nNumOfIters =  this->m_nIterations;

    //Correct for spill-out

    for (int k = 1; k <= nNumOfIters; k++) {

        //Set up image iterator to operate directly on imageEstimate.
        ImageIteratorType it( imageEstimate, imageEstimate->GetLargestPossibleRegion() );
        PixelType voxelVal;

        //Loop over all voxels and remove negative numbers.
        while (! it.IsAtEnd() )
        {
            voxelVal = fmax(it.Get(), 0.0);
            it.Set( voxelVal );
            ++it;
        }

        if ( this->m_bVerbose ) {
            if (k == 1) {
                std::cout << std::endl << "Iteration:  " << std::endl;
            }

            std::cout << k << ":\t";

            i = 0;
            float fNewRegMean = 0;

            labelStatsFilter->SetInput( imageEstimate );
            labelStatsFilter->Update();

            for( typename ValidLabelValuesType::const_iterator vIt=labelStatsFilter->GetValidLabelValues().begin();
                    vIt != labelStatsFilter->GetValidLabelValues().end(); ++vIt)
            {
                if ( labelStatsFilter->HasLabel(*vIt) )
                {
                    LabelPixelType labelValue = *vIt;
                    //std::cout << "Label: " << labelValue << "\t";
                    fNewRegMean = std::max( labelStatsFilter->GetMean( labelValue ), 0.0 );
                    //std::cout << fNewRegMean << " ";
                    vecRegMeansCurrent.put(i, fNewRegMean);
                    i++;
                }
            }

            std::cout << vecRegMeansCurrent << std::endl;
        }

        i=0;

        for( typename ValidLabelValuesType::const_iterator vIt=labelStatsFilter->GetValidLabelValues().begin();
                vIt != labelStatsFilter->GetValidLabelValues().end(); ++vIt)
        {
            if ( labelStatsFilter->HasLabel(*vIt) )
            {

                LabelPixelType labelValue = *vIt;
                binThreshFilter->SetLowerThreshold( labelValue );
                binThreshFilter->SetUpperThreshold( labelValue );

                multiplyFilter->SetInput1( imageEstimate );
                multiplyFilter->SetInput2( binThreshFilter->GetOutput() );

                pBlurFilter->SetInput(multiplyFilter->GetOutput());

                //bkg = bkg + ( cc * (1-mask))

                subtractFilter->SetInput1( 1.0 );
                subtractFilter->SetInput2( binThreshFilter->GetOutput() );

                multiplyFilter2->SetInput1( pBlurFilter->GetOutput() );
                multiplyFilter2->SetInput2( subtractFilter->GetOutput() );
                multiplyFilter2->Update();

                if (i==0){
                    imageBackground = multiplyFilter2->GetOutput();
                }
                else {
                    addFilter->SetInput1(imageBackground);
                    addFilter->SetInput2(multiplyFilter2->GetOutput());
                    addFilter->Update();
                    imageBackground = addFilter->GetOutput();
                }
                imageBackground->DisconnectPipeline();
            }
            i++;
        }

        //- output = ( orig - bkg ) / rec

        subtractFilter->SetInput1( pPET );
        subtractFilter->SetInput2( imageBackground );

        divideFilter->SetInput1( subtractFilter->GetOutput() );
        divideFilter->SetInput2( imageRec );
        divideFilter->Update();

        imageEstimate = divideFilter->GetOutput();
        //imageEstimate->UpdateOutputData();
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
