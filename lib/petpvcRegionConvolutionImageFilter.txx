/*
   petpvcRegionConvolutionPVCImageFilter.txx

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

#ifndef __PETPVCREGIONCONVOLUTIONIMAGEFILTER_TXX
#define __PETPVCREGIONCONVOLUTIONIMAGEFILTER_TXX

#include "petpvcRegionConvolutionImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

using namespace itk;

namespace petpvc
{

template< class TInputImage, class TMaskImage >
RegionConvolutionPVCImageFilter< TInputImage, TMaskImage>
::RegionConvolutionPVCImageFilter()
{
    this->m_bVerbose = false;
}

template< class TInputImage, class TMaskImage >
void RegionConvolutionPVCImageFilter< TInputImage, TMaskImage>
::GenerateData()
{
    typename TInputImage::ConstPointer input = this->GetInput();
    typename TInputImage::Pointer output = this->GetOutput();

    InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));
    MaskImagePointer pMask = dynamic_cast<const TMaskImage*> (ProcessObject::GetInput(1));

    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();    
	typename MultiplyFilterType::Pointer multiplyFilter2 = MultiplyFilterType::New();
    typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();
    typename BlurringFilterType::Pointer blurFilter = BlurringFilterType::New();
	typename BlurringFilterType::Pointer blurFilter2 = BlurringFilterType::New();

	blurFilter->SetVariance( this->GetPSF() );
	blurFilter2->SetVariance( this->GetPSF() );

	//Perform regional convolution
	multiplyFilter->SetInput1( pMask );
	multiplyFilter->SetInput2( pPET );
	blurFilter->SetInput( multiplyFilter->GetOutput() );
	blurFilter2->SetInput( pMask );
	divideFilter->SetInput1( blurFilter->GetOutput() );
	divideFilter->SetInput2( blurFilter2->GetOutput() );
	multiplyFilter2->SetInput1( pMask );
	multiplyFilter2->SetInput2( divideFilter->GetOutput() );
	
	try {
        multiplyFilter2->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "\n[Error]\tCannot performing regional convolution"
                  << std::endl;
    }

	typename TInputImage::Pointer imageOut= multiplyFilter2->GetOutput();
	
    this->AllocateOutputs();

    ImageAlgorithm::Copy( imageOut.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                          output->GetRequestedRegion() );

}

}// end namespace


#endif
