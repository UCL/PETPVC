/*
   petpvcRoussetImageFilter.txx

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

#ifndef __PETPVCROUSSETPVCIMAGEFILTER_HXX
#define __PETPVCROUSSETPVCIMAGEFILTER_HXX
 
#include "petpvcRoussetPVCImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

 
using namespace itk;

namespace petpvc
{
 
template< class TInputImage, class TMaskImage >
void RoussetPVCImageFilter< TInputImage, TMaskImage>
::GenerateData()
{
  typename TInputImage::ConstPointer input = this->GetInput();
  typename TInputImage::Pointer output = this->GetOutput();

	typename GTMImageFilterType::Pointer pGTM = GTMImageFilterType::New();

	InputImagePointer pPET = dynamic_cast<const TInputImage*> (ProcessObject::GetInput(0));
	MaskImagePointer pMask = dynamic_cast<const TMaskImage*> (ProcessObject::GetInput(1));

  pGTM->SetInput( pMask );
  pGTM->SetPSF( this->GetPSF() );
  //Calculate GTM.
  try {
      pGTM->Update();
  } catch (itk::ExceptionObject & err) {
      std::cerr << "[Error]\tCannot calculate GTM"
              << std::endl;
    }

	std::cout << pGTM->GetMatrix() << std::endl;


 
  this->AllocateOutputs();
 
  ImageAlgorithm::Copy(input.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                       output->GetRequestedRegion() );


}
 
}// end namespace
 
 
#endif
