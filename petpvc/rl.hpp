#pragma once

#ifndef _RL_HPP_
#define _RL_HPP_

#include "petpvc.hpp"

namespace  petpvc {

template<typename TInputImage=ImageType3D>
float CalculateLogLikelihood(const typename TInputImage::Pointer a,
                             const typename TInputImage::Pointer b) {

  itk::ImageRegionIterator<TInputImage> it1(a, a->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TInputImage> it2(b, b->GetLargestPossibleRegion());

  float imageLog = 0.0f;

  it1.GoToBegin(); it2.GoToBegin();

  float diff;

  while (!it1.IsAtEnd()){
    if (it2.Get() > 0.0f) {
      diff = it1.Get() * log(it2.Get()) - it2.Get();
      imageLog += diff;
    }

    ++it1; ++it2;
  }

  return imageLog;

}

template<typename TInputImage, typename TBlurFilter>
void RichardsonLucy(const typename TInputImage::Pointer pet,
                const typename TBlurFilter::Pointer blur,
                typename TInputImage::Pointer &output,
                float stoppingCriterion = -3e+06, int niter=10) {

  //TODO: Add 4D processing

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  ImageType3D::Pointer nonNegPET = ImageType3D::New();
  ImageType3D::Pointer currentEstimate = ImageType3D::New();

  //Remove zeros from input.
  RemoveNegativeNumbers(pet,nonNegPET);

  //Set initial estimate to PET without negative numbers.
  Duplicate<ImageType3D>(nonNegPET, currentEstimate);

  int k = 1;
  float logLike = std::numeric_limits<float>::max();

  while ( (k <= niter) && (logLike > stoppingCriterion) ) {

    ImageType3D::Pointer tmp = ImageType3D::New();
    //Smooth current estimate
    ApplySmoothing<ImageType3D,TBlurFilter>(currentEstimate,blur,tmp);
    //Divide non-neg PET by blurred PET
    Divide(nonNegPET,tmp,tmp);
    //Multiply estimate by tmp
    Multiply(currentEstimate,tmp,tmp);
    //Make currentEstimate = tmp
    Duplicate<ImageType3D>(tmp,currentEstimate);
    //Reblur currentEstimate
    ApplySmoothing<ImageType3D,TBlurFilter>(currentEstimate,blur,tmp);
    logLike = CalculateLogLikelihood(nonNegPET,tmp);

    std::cout << k << "\t" << logLike << std::endl;
    k++;
  }

  //Return output
  Duplicate<ImageType3D>(currentEstimate,output);

}

} // end namespace petpvc

#endif