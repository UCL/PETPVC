/*
   FuzzyCorr.cxx

   Author:      Benjamin A. Thomas
 
   Copyright 2013 Institute of Nuclear Medicine, University College London.

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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include "FuzzyCorrFilter.h"

int main(int argc, char *argv[]) {

    //A simple application that applies fuzziness 'correction' and prints
    //the matrix to std out. This could be useful if you simply need the 
    //correction factors to apply as part of some correction technique.

    typedef itk::Image<float, 4> ImageType;
    typedef petpvc::FuzzyCorrFilter<ImageType> FilterType;
    typedef itk::ImageFileReader<ImageType> ReaderType;

    ReaderType::Pointer reader = ReaderType::New();
    FilterType::Pointer filter = FilterType::New();

    if (argc == 2)
        reader->SetFileName(argv[1]);
    else {
        std::cerr << "Usage: pvc_fuzzyCorr <maskfile>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        reader->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;

        return EXIT_FAILURE;
    }

    filter->SetInput(reader->GetOutput());

    try {
        filter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;

        return EXIT_FAILURE;
    }

    vnl_matrix<float> test = filter->GetMatrix();

    test.print(std::cout);


    return EXIT_SUCCESS;
}
