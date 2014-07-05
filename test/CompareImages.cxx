/*=========================================================================
 *  Author:      Kris Thielemans
 *
 *  Copyright University College London
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/*
   This file compares 2 images by taking the absolute difference
   and checking that its maximum is less than a threshold.
*/

#include "itkImageFileReader.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <iostream>


int main( int argc, char *argv[] )
{
    if( argc != 4 ) {
        std::cerr << "Usage: " << argv[0] << " imagefile1 imagefile2 absolute_diff_threshold" << std::endl;
        return EXIT_FAILURE;
    }
    //  We declare the pixel type and dimension of the image to be produced as
    //  output.
    typedef float  PixelType;
    const unsigned int    Dimension = 3;

    typedef itk::Image< PixelType, Dimension >       ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;

    try {
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(argv[1]);
        reader->Update();
        ImageType::Pointer image1 = reader->GetOutput();
        image1->DisconnectPipeline();
        reader->SetFileName(argv[2]);
        reader->Update();
        ImageType::Pointer image2 = reader->GetOutput();
        const float threshold = atof(argv[3]);


        typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType> AbsDiffFilterType;
        AbsDiffFilterType::Pointer absdifffilter = AbsDiffFilterType::New();
        absdifffilter->SetInput1(image1);
        absdifffilter->SetInput2(image2);
        ImageType::Pointer absdiff = absdifffilter->GetOutput();

        absdifffilter->Update();

        typedef itk::MinimumMaximumImageCalculator<ImageType> minmaxCalculatorType;
        minmaxCalculatorType::Pointer minmaxCalculator = minmaxCalculatorType::New();
        minmaxCalculator->SetImage(absdiff);
        minmaxCalculator->ComputeMaximum();

        const float absmax = minmaxCalculator->GetMaximum();
        std::cerr << "max " << absmax << '\n';
        if (absmax > threshold) {
            std::cerr << "Absolute difference too large (" << absmax << ")\n";
            return EXIT_FAILURE;
        }

    } catch( itk::ExceptionObject & excp ) {
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
