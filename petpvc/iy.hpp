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

  ApplySmoothing<TInputImage, TBlurFilter>(pet, blur, output);

  for (int n=1; n < numOfPETVols; n++){
    //For each PET volume

    for (int k=1; i < niter; k++) {

      //Calculate regional means

      //Create synthetic PET
      //Smooth synthetic PET

      //Div s/(s*h)
      //mult pet by div
      //store result for next iter
    }
    //Paste into output volume
  }
  //Return output
}

}; //end namespace petpvc

#endif