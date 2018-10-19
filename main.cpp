#include <iostream>
#include <memory>

#include "petpvc/petpvc.hpp"

using namespace petpvc;

int main(int argc, char *argv[]){

  ImageType3D::Pointer inImage1 = ImageType3D::New();
  MaskImageType3D::Pointer inImage2 = MaskImageType3D::New();
  ImageType3D::Pointer outImage = ImageType3D::New();

  ReadFile<ImageType3D>(argv[1],inImage1);
  ReadFile<MaskImageType3D>(argv[2],inImage2);

  //IterativeYang<ImageType3D, MaskImageType3D>(inImage1, inImage2, outImage);

  BlurringImageFilterType::Pointer gaussian = BlurringImageFilterType::New();
  ConfigureGaussian(gaussian,10,10,10);

  IterativeYang<ImageType3D, MaskImageType3D, BlurringImageFilterType>(inImage1, inImage2, gaussian, outImage);
  //IterativeYang<ImageType3D, MaskImageType4D>(inImage1, inImage2, outImage);
  //IterativeYang<ImageType4D, MaskImageType4D>(inImage1, inImage2, outImage);
  
  WriteFile<ImageType3D>(outImage,"iy-out.nii.gz");


  return EXIT_SUCCESS;
}