#pragma once

#ifndef _GTM_HPP_
#define _GTM_HPP_

#include <iomanip>

#include "petpvc.hpp"

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_matrix_inverse.h>

namespace petpvc {

float CalculateContribution(const ImageType3D::Pointer input,
                            const typename MaskImageType3D::Pointer regionMask) {

  //TODO: Implement for 4D mask.

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
        const typename TBlurFilter::Pointer blur,
        std::vector<float> &outputMeans){

  //TODO: Add 4D processing for PET and masks

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;

  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;

  std::vector<int> labelIndexList;
  GetRegionIndexList(mask,labelIndexList);

  int numOfLabels = labelIndexList.size();

  std::vector<float> regionMeanList;
  GetRegionalMeans(pet,mask,regionMeanList);

  typedef vnl_matrix<float> MatrixType;
  MatrixType weights;
  weights.set_size(numOfLabels, numOfLabels);
  weights.fill(0);

  ImageType3D::Pointer currentVolume = ImageType3D::New();
  MaskImageType3D::Pointer rI = MaskImageType3D::New();
  MaskImageType3D::Pointer rJ = MaskImageType3D::New();
  ImageType3D::Pointer rJsmooth = ImageType3D::New();

  std::vector<float> finalMeans;
  //finalMeans.reserve(numOfLabels);

  std::cout << std::endl;

  for (int n=0; n < numOfPETVols; n++) {
    //For each PET volume
    GetVolume(pet, n, currentVolume);

    for (int i=0; i < numOfLabels; i++){
      //For each region i:

      //Extract region i mask.
      GetRegion(mask,labelIndexList[i],rI);

      //Calculate sum of region I
      float sumOfregionI = CalculateSum<MaskImageType3D>(rI);

      for (int j=0; j < numOfLabels; j++) {
        //For each region j:

        //Extract region j mask
        GetRegion(mask,labelIndexList[j],rJ);
        // Smooth region j
        ApplySmoothing<MaskImageType3D, TBlurFilter>(rJ, blur, rJsmooth);
        //Get contribution of j into i and normalise by sum of region i mask;
        float fracJinI = CalculateContribution(rJsmooth, rI) / sumOfregionI;
        //Put fractional contribution into weights matrix
        weights(i,j) = fracJinI;

      }
    }
    //Print GTM
    std::cout << std::fixed << std::setprecision(4) << weights << std::endl;

    //vnl vector to store the estimated means before GTM correction.
    vnl_vector<float> vnl_regionMeanList;
    vnl_regionMeanList.set_size(numOfLabels);

    for (int n=0; n < numOfLabels; n++){
      vnl_regionMeanList[n] = regionMeanList[n];
    }

    //vnl vector to store the estimated means after GTM correction.
    vnl_vector<float> vnl_regionMeanListUpdated;
    vnl_regionMeanListUpdated.set_size(numOfLabels);

    //Apply matrix inverse to regional mean values.
    vnl_regionMeanListUpdated = vnl_matrix_inverse<float>(weights) * vnl_regionMeanList;

    std::cout << "Corrected means: " << std::fixed << std::setprecision(4) << vnl_regionMeanListUpdated << std::endl;

    for (int n=0; n < numOfLabels; n++){
      finalMeans.push_back(vnl_regionMeanListUpdated[n]);
    }

  }

  outputMeans = finalMeans;
}

} // end namespace petpvc

#endif