#pragma once

#ifndef _IY_HPP
#define _IY_HPP

#include "petpvc.hpp"

namespace petpvc {
template<typename TInputImage, typename TMaskImage, typename TBlurFilter>
void IterativeYang(const typename TInputImage::Pointer pet,
                   const typename TMaskImage::Pointer mask,
                   const typename TBlurFilter::Pointer blur,
                   typename TInputImage::Pointer &output,
                   int niter=10) {

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;

  ImageType3D::Pointer currentVolume = ImageType3D::New();
  ImageType3D::Pointer currentIteration = ImageType3D::New();

  //ApplySmoothing<TInputImage, TBlurFilter>(pet, blur, output);

  for (int n=0; n < numOfPETVols; n++){
    //For each PET volume
    GetVolume(pet,n,currentVolume);

    currentIteration = currentVolume;

    for (int k=1; k < niter; k++) {
      ImageType3D::Pointer tmp = ImageType3D::New();
      //Calculate regional means
      //Create synthetic PET
      //Smooth synthetic PET
      ApplySmoothing<TInputImage, TBlurFilter>(currentIteration, blur, tmp);
      Duplicate<ImageType3D>(tmp, currentIteration);
      //Div s/(s*h)
      //mult pet by div
      //store result for next iter
    }
    //Paste into output volume
  }
  //Return output
  Duplicate<ImageType3D>(currentIteration,output);
}

}; //end namespace petpvc

#endif