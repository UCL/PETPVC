#include <iostream>
#include <memory>

#include "petpvc/petpvc.hpp"

using namespace petpvc;

int main(int argc, char *argv[]){


  //std::shared_ptr<Object> obj1;
  //std::shared_ptr<ImTest3D> obj2(new ImTest3D(argv[1]));

  //auto obj = foo(new ImageObject3D(argv[1]));

  //WriteFile<ImageType3D>(obj2->getImage(),"testWrite.nii.gz");

  //obj1 = obj2;

  //WriteFile<ImageType3D>(obj1->getImage(),"testWriteCast.nii.gz");

  //petpvc::EImageType inputImageType = petpvc::GetImageType(argv[1]);

  //petpvc::ImageType3D::Pointer inImage1 = petpvc::ImageType3D::New();
  //petpvc::ImageType3D::Pointer inImage2 = petpvc::ImageType3D::New();
  //petpvc::ImageType3D::Pointer outImage = petpvc::ImageType3D::New();

  //petpvc::ReadFile<petpvc::ImageType3D>(argv[1],inImage1);
  //petpvc::ReadFile<petpvc::ImageType3D>(argv[2],inImage2);

  //Multiply(inImage1,inImage2,outImage);
  /*
  auto inObj = petpvc::CreateImage(inputImageType, argv[1]);

  for (int i=1; i <= inObj->getNoOfVolumes(); i++) {
    inObj->getVolume(i,inImage2);
    std::stringstream ss;
    ss << "vol_";
    ss << i;
    ss << ".nii.gz";

    petpvc::WriteFile<petpvc::ImageType3D>(inImage2,ss.str());
  }
  */
  /*
  petpvc::ImageType3D::Pointer maskImage = petpvc::ImageType3D::New();

  auto inMaskObj = petpvc::CreateMaskImage(petpvc::EImageType::E4DImage, argv[1]);

  for (int i=1; i <= inMaskObj->getNoOfRegions(); i++) {
    inMaskObj->getRegion(i,maskImage);
    std::stringstream ss;
    ss << "label_";
    ss << i;
    ss << ".nii.gz";

    petpvc::WriteFile<petpvc::ImageType3D>(maskImage,ss.str());
  }*/

  //petpvc::CreateBlankImageFromExample(inImage2,outImage);
  //petpvc::WriteFile<petpvc::MaskImageType3D>(outImage,"blank.nii.gz");

  //std::cout << "Hello, World 2!" << std::endl;

  //petpvc::ImageType4D::Pointer fullImage = petpvc::ImageType4D::New();
  //petpvc::ImageType4D::Pointer blankImage = petpvc::ImageType4D::New();
  //auto inObj = petpvc::CreateImage(petpvc::EImageType::E4DImage, argv[1]);

  //inObj->getImage<petpvc::ImageType4D>(fullImage);

  //std::cout << fullImage;

  //return 0;

  /*
  petpvc::WriteFile<petpvc::ImageType4D>(fullImage,"invert.nii.gz");

  petpvc::CreateBlankImageFromExample<petpvc::ImageType4D>(fullImage, blankImage);

  std::cout << inObj->getNoOfVolumes() << std::endl;

  petpvc::ImageType3D::Pointer vol = petpvc::ImageType3D::New();

  for (int i=1; i < inObj->getNoOfVolumes(); i++){
    int dstPos = inObj->getNoOfVolumes()-i;
    inObj->getVolume(i,vol);
    petpvc::PasteInto(vol,dstPos,blankImage);
  }
   */

  return EXIT_SUCCESS;
}