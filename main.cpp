#include <iostream>
#include <memory>

#include "petpvc/petpvc.hpp"

using namespace petpvc;

int main(int argc, char *argv[]){

  petpvc::EImageType inputImageType = petpvc::GetImageType(argv[1]);

  petpvc::ImageType3D::Pointer inImage1 = petpvc::ImageType3D::New();
  petpvc::ImageType3D::Pointer inImage2 = petpvc::ImageType3D::New();
  petpvc::ImageType3D::Pointer outImage = petpvc::ImageType3D::New();

  //petpvc::ReadFile<petpvc::ImageType3D>(argv[1],inImage1);
  //petpvc::ReadFile<petpvc::ImageType3D>(argv[2],inImage2);

  //Multiply(inImage1,inImage2,outImage);

  auto inObj = petpvc::CreateImage(inputImageType, argv[1]);

  for (int i=1; i <= inObj->getNoOfVolumes(); i++) {
    inObj->getVolume(i,inImage2);
    std::stringstream ss;
    ss << "vol_";
    ss << i;
    ss << ".nii.gz";

    petpvc::WriteFile<petpvc::ImageType3D>(inImage2,ss.str());
  }

  //std::cout << "Hello, World 2!" << std::endl;

  return EXIT_SUCCESS;
}