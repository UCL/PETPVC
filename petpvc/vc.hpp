#pragma once

#ifndef _VC_HPP_
#define _VC_HPP_

#include "petpvc.hpp"

namespace petpvc {

template<typename TInputImage, typename TBlurFilter>
void VanCittert(const typename TInputImage::Pointer pet,
                   const typename TBlurFilter::Pointer blur,
                   typename TInputImage::Pointer &output,
                   float alpha=1.5, float stoppingCriterion = 0.01, int niter=30) {

  //TODO: Add 4D processing

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  ImageType3D::Pointer currentEstimate = ImageType3D::New();
  ImageType3D::Pointer previousEstimate = ImageType3D::New();

  Duplicate<ImageType3D>(pet, currentEstimate);
  //Duplicate(pet, previousEstimate);

  float sumOfPETsq = CalculateSumOfSquares(pet);

  int k = 1;
  float imageDiff = std::numeric_limits<float>::max();

  while ( (k <= niter) && (imageDiff > stoppingCriterion) ) {

    ImageType3D::Pointer tmp = ImageType3D::New();
    //Make previousEstimate = k-1
    Duplicate<ImageType3D>(currentEstimate, previousEstimate);
    //Smooth current estimate
    ApplySmoothing<ImageType3D,TBlurFilter>(currentEstimate,blur,tmp);
    //Subtract smoothed estimate from PET
    Subtract(pet,tmp,tmp);
    //Reblur tmp
    ApplySmoothing<ImageType3D,TBlurFilter>(tmp,blur,tmp);
    //Step by alpha value
    Multiply(tmp,alpha,tmp);
    //Add step to current estimate.
    Add(currentEstimate,tmp,tmp);
    //Enforce non-negativity
    RemoveNegativeNumbers(tmp,tmp);
    //Update estimate
    Duplicate<ImageType3D>(tmp,currentEstimate);
    //Calculate difference
    float sumOfSqDiff = CalculateSumOfSquareDiff(currentEstimate,previousEstimate);
    imageDiff = sqrt( sumOfSqDiff ) / sqrt( sumOfPETsq );
    std::cout << k << "\t" << imageDiff << std::endl;
    k++;
  }

  //Return output
  Duplicate<ImageType3D>(currentEstimate,output);

}

} // end namespace petpvc

#endif