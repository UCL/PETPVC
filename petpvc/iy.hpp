#pragma once

#ifndef _IY_HPP_
#define _IY_HPP_

#include "petpvc.hpp"

namespace petpvc {

void GetSyntheticPETImage(const ImageType3D::Pointer input, const typename MaskImageType3D::Pointer mask,
                                    const std::vector<float> &valList, const std::vector<int> &idxList,
                                    ImageType3D::Pointer &output) {

  ImageType3D::Pointer tmpOutputImage = ImageType3D::New();
  CreateBlankImageFromExample(input,tmpOutputImage);

  for (int n=0; n <idxList.size(); n++){
    ImageType3D::Pointer tmpImage = ImageType3D::New();

    MaskImageType3D::Pointer regionImage = MaskImageType3D::New();
    GetRegion(mask,idxList[n],regionImage);

    typedef itk::CastImageFilter< MaskImageType3D, ImageType3D > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(regionImage);
    castFilter->Update();

    //std::cout <<"\tInserting " << valList[n] << " into " << n << std::endl;
    Multiply(castFilter->GetOutput(),valList[n],tmpImage);

    Add(tmpOutputImage,tmpImage,tmpOutputImage);

  }
  Duplicate<ImageType3D>(tmpOutputImage,output);

}

void GetSyntheticPETImage(const ImageType3D::Pointer input, const typename ImageType4D::Pointer mask,
                                    const std::vector<float> &valList, const std::vector<int> &idxList,
                                    ImageType3D::Pointer &output) {

  //TODO: Implement synthetic PET image generation for 4D.

}

template<typename TInputImage, typename TMaskImage, typename TBlurFilter>
void IterativeYang(const typename TInputImage::Pointer pet,
                   const typename TMaskImage::Pointer mask,
                   const typename TBlurFilter::Pointer blur,
                   typename TInputImage::Pointer &output,
                   int niter=10) {

  //TODO: Add 4D processing

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;

  std::vector<int> labelIndexList;
  GetRegionIndexList(mask,labelIndexList);

  std::vector<float> regionMeanList;

  ImageType3D::Pointer currentVolume = ImageType3D::New();
  ImageType3D::Pointer currentIteration = ImageType3D::New();
  ImageType3D::Pointer synthPET = ImageType3D::New();

  for (int n=0; n < numOfPETVols; n++){
    //For each PET volume
    GetVolume(pet,n,currentVolume);

    currentIteration = currentVolume;

    for (int k=1; k <= niter; k++) {
      std::cout << "Iteration " << k << " : ";
      ImageType3D::Pointer tmp = ImageType3D::New();
      // Calculate regional means
      GetRegionalMeans(currentIteration,mask,regionMeanList);
      // Create synthetic PET
      GetSyntheticPETImage(currentIteration,mask,regionMeanList, labelIndexList,synthPET);
      // Smooth synthetic PET
      ApplySmoothing<ImageType3D, TBlurFilter>(synthPET, blur, tmp);
      // s' = s/[s(x)*h(x)]
      Divide<ImageType3D>(synthPET,tmp,tmp);
      // Multiply pet by s'
      Multiply<ImageType3D>(pet,tmp,tmp);
      // Store result for next iter
      Duplicate<ImageType3D>(tmp, currentIteration);
    }
    //Paste into output volume
  }
  //Return output
  Duplicate<ImageType3D>(currentIteration,output);
}

}; //end namespace petpvc

#endif