#pragma once

#ifndef _LABBE_HPP_
#define _LABBE_HPP_

#include <iomanip>

#include "petpvc.hpp"

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_matrix_inverse.h>

namespace petpvc {

float CalculateContributionSmooth(const ImageType3D::Pointer a,
                            const ImageType3D::Pointer b) {

  //TODO: Implement for 4D mask.

  ImageType3D::Pointer tmpImage = ImageType3D::New();
  Multiply<ImageType3D>(a,b,tmpImage);
  /*
  typedef typename itk::MultiplyImageFilter<ImageType3D> FilterType;
  typename FilterType::Pointer multiply = FilterType::New();
  multiply->SetInput1(input);
  multiply->SetInput2(regionMask);
  multiply->Update();*/


  return CalculateSum<ImageType3D>(tmpImage);

}

template<typename TInputImage, typename TMaskImage, typename TBlurFilter>
void Labbe(const typename TInputImage::Pointer pet,
        const typename TMaskImage::Pointer mask,
        const typename TBlurFilter::Pointer blur,
        std::vector<float> &outputMeans){

  //TODO: Add 4D processing for PET and masks
  //TODO: Add output file format

  const int numOfPETVols = GetNumberOfVolumes<TInputImage>(pet);
  const int numOfMaskVols = GetNumberOfVolumes<TMaskImage>(mask);

#ifndef NDEBUG
  std::cout << "Number of vols. to PV-correct = " << numOfPETVols << std::endl;
  std::cout << "Number of mask vols = " << numOfMaskVols << std::endl;
#endif

  std::vector<int> labelIndexList;
  GetRegionIndexList(mask,labelIndexList);

  int numOfLabels = labelIndexList.size();

  std::vector<float> regionMeanList;
  regionMeanList.reserve(numOfLabels);

  typedef vnl_matrix<float> MatrixType;
  MatrixType weights;
  weights.set_size(numOfLabels, numOfLabels);
  weights.fill(0);

  ImageType3D::Pointer currentVolume = ImageType3D::New();

  MaskImageType3D::Pointer rI = MaskImageType3D::New();
  MaskImageType3D::Pointer rJ = MaskImageType3D::New();

  ImageType3D::Pointer rIsmooth = ImageType3D::New();

  std::vector<float> finalMeans;
  finalMeans.reserve(numOfLabels);

  std::cout << std::endl;

  for (int n=0; n < numOfPETVols; n++) {
    //For each PET volume
    GetVolume(pet, n, currentVolume);

    for (int i=0; i < numOfLabels; i++){
      //For each region i:

      //Extract region i mask.
      GetRegion(mask,labelIndexList[i],rI);

      // Smooth region i
      ApplySmoothing<MaskImageType3D, TBlurFilter>(rI, blur, rIsmooth);

      //Calculate sum of region I
      float sumOfregionI = CalculateSum<ImageType3D>(rIsmooth);

      for (int j=0; j < numOfLabels; j++) {
        //For each region j:
        ImageType3D::Pointer rJsmooth = ImageType3D::New();
        //Extract region j mask
        GetRegion(mask,labelIndexList[j],rJ);
        // Smooth region j
        ApplySmoothing<MaskImageType3D, TBlurFilter>(rJ, blur, rJsmooth);
        //Get contribution of j into i and normalise by sum of region i mask:
        float fSumNeighbour = CalculateContributionSmooth(rJsmooth, rIsmooth);
        float fracJinI = fSumNeighbour / sumOfregionI;
        //Put fractional contribution into weights matrix
        weights(i,j) = fracJinI;

      }

      //Calculate smoothed regional mean.
      ImageType3D::Pointer regPET = ImageType3D::New();
      Multiply(pet,rIsmooth,regPET);
      regionMeanList[i] = CalculateSum(regPET) / CalculateSum(rIsmooth);

    }
    //Print GTM
    std::cout << std::fixed << std::setprecision(6) << weights << std::endl;

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

    std::cout << "Corrected means: " << std::fixed << std::setprecision(6) << vnl_regionMeanListUpdated << std::endl;

    for (int n=0; n < numOfLabels; n++){
      finalMeans.push_back(vnl_regionMeanListUpdated[n]);
    }

  }

  outputMeans = finalMeans;
}

} // end namespace petpvc

#endif