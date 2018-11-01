#include <iostream>
#include <memory>

#include "petpvc/petpvc.hpp"
#include "petpvc/iy.hpp"
#include "petpvc/gtm.hpp"
#include "petpvc/rbv.hpp"
#include "petpvc/vc.hpp"

using namespace petpvc;

int main(int argc, char *argv[]){

  ImageType3D::Pointer inImage1 = ImageType3D::New();
  MaskImageType3D::Pointer inImage2 = MaskImageType3D::New();
  ImageType3D::Pointer outImage = ImageType3D::New();

  ReadFile<ImageType3D>(argv[1],inImage1);
  ReadFile<MaskImageType3D>(argv[2],inImage2);

  BlurringImageFilterType::Pointer gaussian = BlurringImageFilterType::New();
  ConfigureGaussian(gaussian,5,6,7);

  VanCittert<ImageType3D,BlurringImageFilterType>(inImage1,gaussian,outImage,1.5,1e-8,30);
  WriteFile<ImageType3D>(outImage,"vc-3d.nii.gz");

  //GTM<ImageType3D, MaskImageType3D, BlurringImageFilterType>(inImage1, inImage2, gaussian);

  //RBV<ImageType3D, MaskImageType3D, BlurringImageFilterType>(inImage1, inImage2, gaussian, outImage);

  //IterativeYang<ImageType3D, MaskImageType3D, BlurringImageFilterType>(inImage1, inImage2, gaussian, outImage,50);
  //IterativeYang<ImageType3D, MaskImageType4D>(inImage1, inImage2, outImage);
  //IterativeYang<ImageType4D, MaskImageType4D>(inImage1, inImage2, outImage);
  
  //WriteFile<ImageType3D>(outImage,"rbv-3d.nii.gz");

  /*
  std::vector<int> labelIdx;
  GetRegionIndexList(inImage2,labelIdx);

  std::vector<float> valList;
  GetRegionalMeans(inImage1,inImage2,valList);

  std::cout << "valList = ";
  for (auto x : valList){
    std::cout << x << " ";
  }
  std::cout << std::endl;

  std::cout << "labelIdx = ";
  for (auto x : labelIdx){
    std::cout << x << " ";
  }
  std::cout << std::endl;

  GetSyntheticPETImage(inImage1,inImage2,valList,labelIdx,outImage);
  WriteFile<ImageType3D>(outImage,"iy-synth.nii.gz");*/





  return EXIT_SUCCESS;
}