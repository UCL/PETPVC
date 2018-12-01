#pragma once

#ifndef _MG_HPP_
#define _MG_HPP_

#include "petpvc.hpp"

namespace petpvc {

template<typename TInputImage, typename TMaskImage, typename TBlurFilter>
bool MullerGartner(const typename TInputImage::Pointer pet,
                   const typename TMaskImage::Pointer mask,
                   const typename TBlurFilter::Pointer blur,
                   typename TInputImage::Pointer &output) {

  //TODO: Add 4D processing

  std::cout << "Muller-Gartner" << std::endl;

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;

  std::vector<int> labelIndexList;
  GetRegionIndexList(mask, labelIndexList);


  if (labelIndexList.size() < 3){
    std::cerr << "A minimum of 2 labels (GM + WM) is required for correction!";
    return false;
  }

  typename TInputImage::Pointer currPETImage = TInputImage::New();
  Duplicate<TInputImage>(pet,currPETImage);

  //Get GM
  MaskImageType3D::Pointer gmImage = MaskImageType3D::New();
  GetRegion(mask,labelIndexList[0],gmImage);
  //Get WM
  MaskImageType3D::Pointer wmImage = MaskImageType3D::New();
  GetRegion(mask,labelIndexList[1],wmImage);
  //(Get CSF)
  if (labelIndexList.size()-1 == 3){
    std::cout << "Has CSF image" << std::endl;
    //currPETImage
  }

  //Erode WM
  //PET * Eroded WM
  //Calc mean or replace with given WM value
  //WM * WM value = vWM
  //Smooth vWM
  //Subtract smooth vWM from PET - vGM
  //Clip vGM to GM only (Check order)

  //Smooth GM
  //Divide vGM by sGM
  //Clean-up infs / nans

  return true;
}

} //end namespace petpvc

#endif