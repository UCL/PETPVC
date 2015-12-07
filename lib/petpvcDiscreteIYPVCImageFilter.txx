/*
   petpvcDiscreteIYPVCImageFilter.txx

   Author:      Benjamin A. Thomas

   Copyright 2014-2015 Institute of Nuclear Medicine, University College London.
   Copyright 2014-2015 Clinical Imaging Research Centre, A*STAR-NUS.

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

#ifndef __PETPVCDiscreteIYPVCIMAGEFILTER_TXX
#define __PETPVCDiscreteIYPVCIMAGEFILTER_TXX

#include "petpvcDiscreteIYPVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
DiscreteIYPVCImageFilter< TInputImage, TMaskImage>
::DiscreteIYPVCImageFilter()
{
    this->m_nIterations = 10;
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void DiscreteIYPVCImageFilter< TInputImage, TMaskImage>
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

    typename TInputImage::Pointer imageExtractedRegion;// = InputImagePointer::New();

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

    //For applying the Yang correction step:

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

        labelStatsFilter->SetInput( imageEstimate );
        labelStatsFilter->Update();

        if ( this->m_bVerbose ) {
            if (k == 1) {
                std::cout << std::endl << "Iteration:  " << std::endl;
            }

            std::cout << k << ":\t"; // << std::endl;
            //std::flush(std::cout);
        }

        int i = 0;
        float fNewRegMean = 0;

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

        vecRegMeansUpdated = vecRegMeansCurrent;

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

        i=0;

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

                //Multiply current image estimate by region mask. To clip PET values
                //to mask.
                multiplyFilter->SetInput1( vecRegMeansUpdated.get(i) );
                multiplyFilter->SetInput2( imageExtractedRegion );
                multiplyFilter->Update();

                //If this is the first region, create imageYang,
                //else add the current region to the previous contents of imageYang.
                if (i == 0) {
                    imageYang = multiplyFilter->GetOutput();
                    imageYang->DisconnectPipeline();
                } else {
                    addFilter->SetInput1(imageYang);
                    addFilter->SetInput2(multiplyFilter->GetOutput());
                    addFilter->Update();

                    imageYang = addFilter->GetOutput();
                }

                i++;
            }
        }        

        //Takes the original PET data and the pseudo PET image, calculates the
        //correction factors  and returns the PV-corrected PET image.

        pBlurFilter->SetInput(imageYang);
        pBlurFilter->SetVariance( this->GetPSF() );

        //Take ratio of pseudo PET and smoothed pseudo PET. These are the correction
        //factors.
        divideFilter->SetInput1( imageYang );
        divideFilter->SetInput2( pBlurFilter->GetOutput() );

        //Multiply original PET by correction factors.
        multiplyFilter2->SetInput1( pPET );
        multiplyFilter2->SetInput2( divideFilter->GetOutput() );
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
