#include <iostream>
#include <memory>

#include "itkAddImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSubtractImageFilter.h"

namespace petpvc {

typedef itk::Image<float, 3> ImageType3D;
typedef itk::Image<float, 4> ImageType4D;

typedef itk::Image<short, 3> MaskImageType3D;

enum class EImageType { E3DImage, E4DImage, EUnknown };

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

template<typename TImageType=ImageType3D, typename TOutputImage=TImageType>
void CreateBlankImageFromExample(const typename TImageType::Pointer input, typename TOutputImage::Pointer &output)
{
  typedef itk::ImageDuplicator< TImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(input);
  duplicator->Update();

  typedef itk::CastImageFilter< TImageType, TOutputImage > CastImageFilterType;
  typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
  castFilter->SetInput( duplicator->GetOutput() );
  castFilter->Update();
  typename TOutputImage::Pointer clonedImage = castFilter->GetOutput();

  typedef typename TOutputImage::PixelType PixelType;
  clonedImage->FillBuffer( itk::NumericTraits< PixelType >::Zero );

  output->Graft(clonedImage);

}

template<typename TImageType=ImageType4D, typename TOutputImage=ImageType3D>
void Extract4DVolTo3D(const typename TImageType::Pointer input,
                     const int n,
                     typename TOutputImage::Pointer &output) {

  std::cout << "4D image extraction: " << std::endl;

  typename TImageType::IndexType desiredStart;
  desiredStart.Fill(0);

  typename TImageType::SizeType desiredSize =
      input->GetLargestPossibleRegion().GetSize();

  if (n - 1 < 0) {
    std::cerr << "Requested region " << n << " out of range!" << std::endl;
    throw false;
  }

  desiredStart[3] = n - 1;
  desiredSize[3] = 0;

  typedef itk::ExtractImageFilter<TImageType, TOutputImage> ExtractFilterType;
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();

  extractFilter->SetExtractionRegion(
      typename ImageType4D::RegionType(desiredStart, desiredSize));
  extractFilter->SetInput(input);
  extractFilter->SetDirectionCollapseToIdentity();
  extractFilter->Update();

  output->Graft(extractFilter->GetOutput());

}

class PETPVCImageObject {

public:
  PETPVCImageObject(){};
  explicit PETPVCImageObject(const std::string &filename);
  virtual void getVolume( const int n, ImageType3D::Pointer &vol )=0;
  virtual int getNoOfVolumes()=0;
  virtual ~PETPVCImageObject(){};

protected:
  std::string _inFilename;
};

PETPVCImageObject::PETPVCImageObject(const std::string &filename){
  _inFilename = filename;
  std::cout << "Base ctor" << std::endl;
};

class PETPVCMaskObject : public PETPVCImageObject {
public:
  virtual void getRegion( const int n, ImageType3D::Pointer &reg )=0;
  virtual void getRegion( const int n, MaskImageType3D::Pointer &reg )=0;
  virtual int getNoOfRegions(){ return _numOfRegions; };

protected:
  int _numOfRegions = 0;
};

class ImageIn4D : public PETPVCImageObject {
public:
  ImageIn4D(const std::string &filename);

  void getVolume( const int n, ImageType3D::Pointer &vol );
  int getNoOfVolumes(){ return _internalImage->GetLargestPossibleRegion().GetSize(3); };
protected:
  ImageType4D::Pointer _internalImage;
};

class ImageIn3D : public PETPVCImageObject {
public:
  ImageIn3D(){};
  ImageIn3D(const std::string &filename);

  void getVolume( const int n, ImageType3D::Pointer &vol );
  int getNoOfVolumes(){ return 1; };

protected:
  ImageType3D::Pointer _internalImage;
};

class MaskIn3D : public PETPVCMaskObject {
public:
  MaskIn3D(const std::string &filename);

  //For getting information about the mask labels
  typedef itk::LabelStatisticsImageFilter<ImageType3D, MaskImageType3D> LabelStatisticsFilterType;
  typedef typename LabelStatisticsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  typedef typename LabelStatisticsFilterType::LabelPixelType LabelPixelType;

  typedef itk::BinaryThresholdImageFilter<MaskImageType3D, MaskImageType3D> BinaryThresholdImageFilterType;

  void getVolume( const int n, ImageType3D::Pointer &vol );
  int getNoOfVolumes(){ return 1; };

  void getRegion( const int n, ImageType3D::Pointer &reg );
  void getRegion( const int n, MaskImageType3D::Pointer &reg );

protected:
  MaskImageType3D::Pointer _internalImage;
  std::vector<int> _LabelIndices;

};

class MaskIn4D : public PETPVCMaskObject {
public:
  MaskIn4D(const std::string &filename);

  void getVolume( const int n, ImageType3D::Pointer &vol );
  int getNoOfVolumes(){ return _internalImage->GetLargestPossibleRegion().GetSize(3); };
  void getRegion( const int n, ImageType3D::Pointer &reg );
  void getRegion( const int n, MaskImageType3D::Pointer &reg );
protected:
  ImageType4D::Pointer _internalImage;
};

ImageIn3D::ImageIn3D(const std::string &filename){
  _inFilename = filename;
  _internalImage = ImageType3D::New();
  ReadFile<ImageType3D>( _inFilename, _internalImage );
}

ImageIn4D::ImageIn4D(const std::string &filename){
  _inFilename = filename;
  _internalImage = ImageType4D::New();
  ReadFile<ImageType4D>( _inFilename, _internalImage );
}

MaskIn3D::MaskIn3D(const std::string &filename){

  _inFilename = filename;
  _internalImage = MaskImageType3D::New();
  ReadFile<MaskImageType3D>( _inFilename, _internalImage );

  //Get labels
  petpvc::ImageType3D::Pointer tmpImage = petpvc::ImageType3D::New();
  CreateBlankImageFromExample<MaskImageType3D, ImageType3D>(_internalImage, tmpImage);

  std::cout << "Created tmpImage" << std::endl;

  typename LabelStatisticsFilterType::Pointer labelStatsFilter = LabelStatisticsFilterType::New();
  labelStatsFilter->SetInput(tmpImage);
  labelStatsFilter->SetLabelInput(_internalImage);
  labelStatsFilter->Update();

  _numOfRegions = labelStatsFilter->GetNumberOfLabels()-1;

  std::cout << "Number of labels found: " << _numOfRegions << std::endl;

  //TODO: Die if no labels found.

  //Put labels into array/vector
  _LabelIndices.reserve(_numOfRegions+1);

  std::cout << "Found Label(s): ";
  for( typename ValidLabelValuesType::const_iterator vIt=labelStatsFilter->GetValidLabelValues().begin();
       vIt != labelStatsFilter->GetValidLabelValues().end(); ++vIt) {
    if (labelStatsFilter->HasLabel(*vIt)) {
      LabelPixelType labelValue = *vIt;
      std::cout << labelValue << " ";
      _LabelIndices.push_back(labelValue);
    }
  }

  std::cout << std::endl;
  std::cout << "Label indices: " ;

  for (auto x : _LabelIndices) {
    std::cout <<  x << " ";
  }
  std::cout << std::endl;

}

MaskIn4D::MaskIn4D(const std::string &filename) {

  _inFilename = filename;
  _internalImage = ImageType4D::New();
  ReadFile<ImageType4D>(_inFilename, _internalImage);

  _numOfRegions = _internalImage->GetLargestPossibleRegion().GetSize(3);

}

void ImageIn3D::getVolume( const int n, ImageType3D::Pointer &vol ) {
  std::cout << "3D image: " << _inFilename << std::endl;

  if ( n != 1) {
    std::cerr << "Requested region " << n << " out of range!" << std::endl;
    throw false;
  }
  vol->Graft(_internalImage);
}

void MaskIn3D::getVolume( const int n, ImageType3D::Pointer &vol ) {
  std::cout << "3D mask image: " << _inFilename << std::endl;

  if ( n != 1) {
    std::cerr << "Requested region " << n << " out of range!" << std::endl;
    throw false;
  }

  typedef itk::CastImageFilter< MaskImageType3D, ImageType3D > CastImageFilterType;
  typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
  castFilter->SetInput( _internalImage );
  castFilter->Update();

  vol->Graft(castFilter->GetOutput());
}


void ImageIn4D::getVolume( const int n, ImageType3D::Pointer &vol ) {
  std::cout << "4D image: " << _inFilename << std::endl;

  Extract4DVolTo3D(_internalImage, n, vol);

};

void MaskIn3D::getRegion( const int n, MaskImageType3D::Pointer &reg ){

  typename BinaryThresholdImageFilterType::Pointer binThreshFilter = BinaryThresholdImageFilterType::New();

  int targetIdx = _LabelIndices[n];

  std::cout << "Extracting region " << targetIdx << std::endl;

  binThreshFilter->SetInput( _internalImage );
  binThreshFilter->SetInsideValue(1);
  binThreshFilter->SetOutsideValue(0);

  binThreshFilter->SetLowerThreshold( targetIdx );
  binThreshFilter->SetUpperThreshold( targetIdx );
  binThreshFilter->Update();

  reg->Graft(binThreshFilter->GetOutput());

}

void MaskIn3D::getRegion( const int n, ImageType3D::Pointer &reg ) {

  MaskImageType3D::Pointer tmpRegion = MaskImageType3D::New();

  this->getRegion(n, tmpRegion);

  typedef itk::CastImageFilter< MaskImageType3D, ImageType3D > CastImageFilterType;
  typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
  castFilter->SetInput( tmpRegion );
  castFilter->Update();

  reg->Graft(castFilter->GetOutput());
}
void MaskIn4D::getVolume( const int n, ImageType3D::Pointer &vol ) {
  Extract4DVolTo3D(_internalImage,n,vol);
}

void MaskIn4D::getRegion(const int n, ImageType3D::Pointer &reg){
  Extract4DVolTo3D(_internalImage,n,reg);
}

void MaskIn4D::getRegion(const int n, MaskImageType3D::Pointer &reg){

  ImageType3D::Pointer tmpRegion = ImageType3D::New();
  this->getRegion(n,tmpRegion);

  typedef itk::CastImageFilter< ImageType3D, MaskImageType3D > CastImageFilterType;
  typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
  castFilter->SetInput( tmpRegion );
  castFilter->Update();

  reg->Graft(castFilter->GetOutput());

}

std::unique_ptr<PETPVCImageObject> CreateImage(EImageType e, const std::string &filename){
  if (e == EImageType::E3DImage)
    return std::unique_ptr<PETPVCImageObject>(new ImageIn3D(filename));
  if (e == EImageType::E4DImage)
    return std::unique_ptr<PETPVCImageObject>(new ImageIn4D(filename));

  return nullptr;
}

std::unique_ptr<PETPVCMaskObject> CreateMaskImage(EImageType e, const std::string &filename){
  if (e == EImageType::E3DImage)
    return std::unique_ptr<PETPVCMaskObject>(new MaskIn3D(filename));
  if (e == EImageType::E4DImage)
    return std::unique_ptr<PETPVCMaskObject>(new MaskIn4D(filename));

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