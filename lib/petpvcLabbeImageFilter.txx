/*
   petpvcLabbeImageFilter.txx

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

#ifndef __PETPVCLABBEIMAGEFILTER_TXX
#define __PETPVCLABBEIMAGEFILTER_TXX

#include "petpvcLabbeImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <itkStatisticsImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkImageDuplicator.h>
#include "vnl/vnl_matrix.h"

using namespace itk;

namespace petpvc
{

template<class TImage>
LabbeImageFilter<TImage>::LabbeImageFilter()
{
    //Constructor. Just initialises the matrix and vector that will
    //contain the results of the correction that the filter implements.
    this->matCorrFactors = new MatrixType;
    this->vecSumOfRegions = new VectorType;

    //this->vecVariance = new ITKVectorType;

}

template<class TImage>
void LabbeImageFilter<TImage>::GenerateData()
{

    //Get pointers to input and output.
    typename TImage::ConstPointer input = this->GetInput();

    //Get region size.
    typename TImage::SizeType imageSize =
        input->GetLargestPossibleRegion().GetSize();

    int nClasses = 0;

    //Set size of output matrix and vector.
    //TODO: Should throw an exception here if imageSize.Dimension != 4.
    if (imageSize.Dimension == 4) {
        nClasses = imageSize[3];
        matCorrFactors->set_size(nClasses, nClasses);
        matCorrFactors->fill(0);

        vecSumOfRegions->set_size(nClasses);
        vecSumOfRegions->fill(0);

    }

    typedef itk::Image<float, 3> MaskImageType;

    typedef itk::StatisticsImageFilter<MaskImageType> StatisticsFilterType;
    typedef itk::ExtractImageFilter<TImage, MaskImageType> ExtractFilterType;
    typedef itk::MultiplyImageFilter<MaskImageType, MaskImageType> MultiplyFilterType;
    typedef itk::DiscreteGaussianImageFilter<MaskImageType, MaskImageType> BlurringFilterType;

    typedef itk::ImageDuplicator<MaskImageType> DuplicatorType;

    typename TImage::IndexType desiredStart;
    desiredStart.Fill(0);

    typename TImage::SizeType desiredSize =
        input->GetLargestPossibleRegion().GetSize();

    typename TImage::RegionType desiredRegion;

    StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();

    typename ExtractFilterType::Pointer extractTargetFilter =
        ExtractFilterType::New();
    typename ExtractFilterType::Pointer extractNeighbourFilter =
        ExtractFilterType::New();
    typename MultiplyFilterType::Pointer multiplyFilter =
        MultiplyFilterType::New();
    typename BlurringFilterType::Pointer blurringFilter =
        BlurringFilterType::New();

    blurringFilter->SetVariance((this->GetPSF()));

	typename BlurringFilterType::Pointer blurringFilter2 =
        BlurringFilterType::New();

    blurringFilter2->SetVariance((this->GetPSF()));

    float fSumTarget;
    float fSumNeighbour;

    for (int i = 1; i <= nClasses; i++) {

        fSumTarget = 0.0;
        desiredStart[3] = i - 1;
        desiredSize[3] = 0;

        //Extract 3D brain mask volume i from 4D image.
        extractTargetFilter->SetExtractionRegion(
            typename TImage::RegionType(desiredStart, desiredSize));
        extractTargetFilter->SetInput(input);
        extractTargetFilter->SetDirectionCollapseToIdentity(); // This is required.
        extractTargetFilter->Update();

        blurringFilter->SetInput(extractTargetFilter->GetOutput());

        statsFilter->SetInput(blurringFilter->GetOutput());
        statsFilter->Update();

        //Calculate the sum of non-zero voxels.
        fSumTarget = statsFilter->GetSum();

        vecSumOfRegions->put(i - 1, fSumTarget);

        for (int j = 1; j <= nClasses; j++) {
            fSumNeighbour = 0.0;
            desiredStart[3] = j - 1;
            desiredSize[3] = 0;

            //Extract 3D brain mask volume j from 4D image.
            extractNeighbourFilter->SetExtractionRegion(
                typename TImage::RegionType(desiredStart, desiredSize));
            extractNeighbourFilter->SetInput(input);
            extractNeighbourFilter->SetDirectionCollapseToIdentity(); // This is required.
            extractNeighbourFilter->Update();

			blurringFilter2->SetInput(extractNeighbourFilter->GetOutput());

            //Multiply i by j.
            multiplyFilter->SetInput1(blurringFilter->GetOutput());
            multiplyFilter->SetInput2(blurringFilter2->GetOutput());

            statsFilter->SetInput(multiplyFilter->GetOutput());
            statsFilter->Update();

            //Calculate the sum of the remaining voxels.
            fSumNeighbour = statsFilter->GetSum();

            //Fill location in matrix with sum of remaining voxels
            //normalised by size of i.
            matCorrFactors->put(i - 1, j - 1, fSumNeighbour / fSumTarget);

            //std::cout << "Target: " << fSumTarget << " Neighbour: " << fSumNeighbour << std::endl;
        }

    }
}

template<class TImage>
LabbeImageFilter<TImage>::~LabbeImageFilter()
{
    //Destructor.
    delete this->matCorrFactors;
    delete this->vecSumOfRegions;

}

} // end namespace

#endif
