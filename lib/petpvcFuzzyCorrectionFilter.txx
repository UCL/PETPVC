/*
   petpvcFuzzyCorrectionFilter.txx

   Authors:     Benjamin A. Thomas
                Kris Thielemans (minor modifications)

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

#ifndef __PETPVCFUZZYCORRECTIONFILTER_TXX
#define __PETPVCFUZZYCORRECTIONFILTER_TXX

#include "petpvcFuzzyCorrectionFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <itkStatisticsImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include "vnl/vnl_matrix.h"

using namespace itk;

namespace petpvc
{

template<class TImage>
FuzzyCorrectionFilter<TImage>::FuzzyCorrectionFilter()
{
    //Constructor. Just initialises the matrix and vector that will
    //contain the results of the correction that the filter implements.
    this->matFuzz = new MatrixType;
    this->vecSumOfRegions = new VectorType;
}

template<class TImage>
void FuzzyCorrectionFilter<TImage>::GenerateData()
{

    //Get pointers to input and output.
    typename TImage::ConstPointer input = this->GetInput();
    typename TImage::Pointer output = this->GetOutput();

    this->AllocateOutputs();

    //Get region size.
    typename TImage::SizeType imageSize =
        input->GetLargestPossibleRegion().GetSize();

    int nClasses = 0;

    //Set size of output matrix and vector.
    //TODO: Should throw an exception here if imageSize.Dimension != 4.
    if (imageSize.Dimension == 4) {
        nClasses = imageSize[3];
        matFuzz->set_size(nClasses, nClasses);
        matFuzz->fill(0);

        vecSumOfRegions->set_size(nClasses);
        vecSumOfRegions->fill(0);

    }

    typedef itk::Image<float, 3> MaskImageType;

    typedef itk::StatisticsImageFilter<MaskImageType> StatisticsFilterType;
    typedef itk::ExtractImageFilter<TImage, MaskImageType> ExtractFilterType;
    typedef itk::MultiplyImageFilter<MaskImageType, MaskImageType> MultiplyFilterType;

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

        statsFilter->SetInput(extractTargetFilter->GetOutput());
        statsFilter->Update();

        //Calculate the sum of non-zero voxels.
        fSumTarget = statsFilter->GetSum();
        if (fSumTarget == 0) {
          itkExceptionMacro("Region " << i << " has zero sum, i.e. no voxels in mask. Remove this region.");
        }
        if (statsFilter->GetMinimum() < 0) {
          itkExceptionMacro("Region " << i << "contains negative voxels in mask. Remove this region.");
        }
        std::cerr << "sum in region " << i << " = " << fSumTarget << std::endl;

        vecSumOfRegions->put(i - 1, fSumTarget);

        //Set region i as first input.
        multiplyFilter->SetInput1(extractTargetFilter->GetOutput());

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

            //Multiply i by j.
            multiplyFilter->SetInput2(extractNeighbourFilter->GetOutput());

            statsFilter->SetInput(multiplyFilter->GetOutput());
            statsFilter->Update();

            //Calculate the sum of the remaining voxels.
            fSumNeighbour = statsFilter->GetSum();

            //Fill location in matrix with sum of remaining voxels
            //normalised by size of i.
            matFuzz->put(i - 1, j - 1, fSumNeighbour / fSumTarget);

        }

    }
}

template<class TImage>
FuzzyCorrectionFilter<TImage>::~FuzzyCorrectionFilter(void)
{
    //Destructor.
    delete this->matFuzz;
    delete this->vecSumOfRegions;
}
} // end namespace

#endif
