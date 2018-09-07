#include <iostream>
#include <memory>

#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSubtractImageFilter.h"

namespace petpvc {

typedef itk::Image<float, 3> ImageType3D;
typedef itk::Image<float, 4> ImageType4D;

enum class EImageType { E3DImage, E4DImage, EUnknown };

/*
class PETPVCImageObject {

public:
  std::string getFilename(){ return _filename; };

protected:
  std::string _filename;
  virtual void getVolume( const int n, ImageType3D::Pointer &vol )=0;
  PETPVCImageObject *obj;

};*/

class PETPVCImageObject {

public:
  PETPVCImageObject(std::string filename);
  virtual void getVolume( const int n, ImageType3D::Pointer &vol )=0;
  virtual ~PETPVCImageObject(){};

protected:
  std::string _inFilename;
};

PETPVCImageObject::PETPVCImageObject(std::string filename){
  _inFilename = filename;
  std::cout << "Base ctor" << std::endl;
};

class ImageIn4D : public PETPVCImageObject {
  using PETPVCImageObject::PETPVCImageObject;
  void getVolume( const int n, ImageType3D::Pointer &vol ) { std::cout << "4D image: " << _inFilename << std::endl; };

protected:
  ImageType4D::Pointer _internalImage;
};

class ImageIn3D : public PETPVCImageObject {
  using PETPVCImageObject::PETPVCImageObject;
  void getVolume( const int n, ImageType3D::Pointer &vol ) { std::cout << "3D image: " << _inFilename << std::endl; };

protected:
  ImageType3D::Pointer _internalImage;
};

std::unique_ptr<PETPVCImageObject> CreateImage(EImageType e, const std::string &filename){
  if (e == EImageType::E3DImage)
    return std::unique_ptr<PETPVCImageObject>(new ImageIn3D(filename));
  else if (e == EImageType::E4DImage)
    return std::unique_ptr<PETPVCImageObject>(new ImageIn4D(filename));
  else
    return nullptr;
}

static bool GetImageIO(const std::string &filename, itk::ImageIOBase::Pointer &imageIO) {

  imageIO = itk::ImageIOFactory::CreateImageIO(
      filename.c_str(), itk::ImageIOFactory::ReadMode);

  if (!imageIO) {
    std::cerr << "Could not CreateImageIO for: " << filename << std::endl;
    return false;
  }

  imageIO->SetFileName(filename);
  imageIO->ReadImageInformation();

  return true;
}

static EImageType GetImageType(const std::string &filename){

  itk::ImageIOBase::Pointer imageIO;
  petpvc::GetImageIO(filename, imageIO);

  const size_t numOfDims = imageIO->GetNumberOfDimensions();

  std::cout << "numDimensions: " << numOfDims << std::endl;

  switch ( numOfDims ){
    case 3 : return EImageType::E3DImage; break;
    case 4 : return EImageType::E4DImage; break;
  }

  return EImageType::EUnknown;
}

template<typename TImageType>
void ReadFile(const std::string &filename, typename TImageType::Pointer image)
{
  typedef itk::ImageFileReader<TImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(filename);
  reader->Update();

  image->Graft(reader->GetOutput());
}

template<typename TImageType>
void WriteFile(typename TImageType::Pointer image, const std::string &filename)
{
  typedef typename itk::ImageFileWriter<TImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetFileName(filename);
  writer->Update();
}

template<typename TImageType=ImageType3D>
void Add(const typename TImageType::Pointer a,
         const typename TImageType::Pointer b,
         typename TImageType::Pointer outputImage)
{
  typedef typename itk::AddImageFilter<TImageType> FilterType;
  typename FilterType::Pointer add = FilterType::New();
  add->SetInput1(a);
  add->SetInput2(b);
  add->Update();

  outputImage->Graft(add->GetOutput());
}

template<typename TImageType=ImageType3D>
void Subtract(const typename TImageType::Pointer a,
              const typename TImageType::Pointer b,
              typename TImageType::Pointer outputImage)
{
  typedef typename itk::SubtractImageFilter<TImageType> FilterType;
  typename FilterType::Pointer sub = FilterType::New();
  sub->SetInput1(a);
  sub->SetInput2(b);
  sub->Update();

  outputImage->Graft(sub->GetOutput());
}

template<typename TImageType=ImageType3D>
void Multiply(const typename TImageType::Pointer a,
              const typename TImageType::Pointer b,
              typename TImageType::Pointer outputImage)
{
  typedef typename itk::MultiplyImageFilter<TImageType> FilterType;
  typename FilterType::Pointer multiply = FilterType::New();
  multiply->SetInput1(a);
  multiply->SetInput2(b);
  multiply->Update();

  outputImage->Graft(multiply->GetOutput());
}

template<typename TImageType=ImageType3D>
void Divide(const typename TImageType::Pointer a,
            const typename TImageType::Pointer b,
            typename TImageType::Pointer outputImage)
{
  typedef typename itk::MultiplyImageFilter<TImageType> FilterType;
  typename FilterType::Pointer divide = FilterType::New();
  divide->SetInput1(a);
  divide->SetInput2(b);
  divide->Update();

  outputImage->Graft(divide->GetOutput());
}

}