#pragma once

#ifndef _RBV_HPP_
#define _RBV_HPP_

#include "petpvc.hpp"
#include "gtm.hpp"
#include "iy.hpp"

namespace petpvc {


template<typename TInputImage, typename TMaskImage, typename TBlurFilter>
void RBV(const typename TInputImage::Pointer pet,
                   const typename TMaskImage::Pointer mask,
                   const typename TBlurFilter::Pointer blur,
                   typename TInputImage::Pointer &output) {

  //TODO: Add 4D processing

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;

  std::vector<int> labelIndexList;
  GetRegionIndexList(mask, labelIndexList);

  std::vector<float> gtmMeans;
  GTM<TInputImage, TMaskImage, TBlurFilter>(pet,mask,blur,gtmMeans);

  std::cout << "GTM corrected means: ";

  for (int n=0; n < gtmMeans.size(); n++){
   std::cout << gtmMeans[n] << " ";
  }
  std::cout << std::endl;

  ImageType3D::Pointer synthPET = ImageType3D::New();
  GetSyntheticPETImage(pet,mask,gtmMeans,labelIndexList,synthPET);

  ImageType3D::Pointer synthPETsmooth = ImageType3D::New();
  ApplySmoothing<ImageType3D, TBlurFilter>(synthPET,blur,synthPETsmooth);

  ImageType3D::Pointer ratio = ImageType3D::New();
  Divide(synthPET,synthPETsmooth,ratio);

  Multiply(pet,ratio,output);


}
} //end namespace petpvc

#endif