#pragma once

#ifndef _GTM_HPP_
#define _GTM_HPP_

#include "petpvc.hpp"

namespace petpvc {

template<typename TInputImage=ImageType3D>
float CalculateSum(const typename TInputImage::Pointer input) {

  itk::ImageRegionIterator<TInputImage> it(input, input->GetLargestPossibleRegion());

  float imageSum = 0.0f;
  for (it.Begin(); !it.IsAtEnd(); ++it) {
    imageSum += it.Get();
  }

  return imageSum;

}
float CalculateContribution(const ImageType3D::Pointer input,
                            const typename MaskImageType3D::Pointer regionMask) {

  typedef itk::CastImageFilter< MaskImageType3D, ImageType3D > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(regionMask);
  castFilter->Update();

  ImageType3D::Pointer tmpImage = ImageType3D::New();
  Multiply(input,castFilter->GetOutput(),tmpImage);

  return CalculateSum(tmpImage);

}

template<typename TInputImage, typename TMaskImage, typename TBlurFilter>
void GTM(const typename TInputImage::Pointer pet,
        const typename TMaskImage::Pointer mask,
        const typename TBlurFilter::Pointer blur){

  //TODO: Add 4D processing
  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;

  std::vector<int> labelIndexList;
  GetRegionIndexList(mask,labelIndexList);

  int numOfLabels = labelIndexList.size();

  std::vector<float> regionMeanList;

  ImageType3D::Pointer currentVolume = ImageType3D::New();
  ImageType3D::Pointer currentIteration = ImageType3D::New();
  ImageType3D::Pointer synthPET = ImageType3D::New();

  std::cout << std::endl;

  for (int n=0; n < numOfPETVols; n++) {
    //For each PET volume
    GetVolume(pet, n, currentVolume);

    currentIteration = currentVolume;

    for (int i=0; i < numOfLabels; i++){

      MaskImageType3D::Pointer rI = MaskImageType3D::New();
      GetRegion(mask,labelIndexList[i],rI);

      float sumOfregionI = CalculateSum<MaskImageType3D>(rI);

      for (int j=0; j < numOfLabels; j++) {
        ImageType3D::Pointer rJsmooth = ImageType3D::New();

        MaskImageType3D::Pointer rJ = MaskImageType3D::New();
        GetRegion(mask,labelIndexList[j],rJ);

        // Smooth region j
        ApplySmoothing<MaskImageType3D, TBlurFilter>(rJ, blur, rJsmooth);

        float fracJinI = CalculateContribution(rJsmooth, rI) / sumOfregionI;

        std::cout << fracJinI << " ";

      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;


}


} // end namespace petpvc

#endif