#pragma once

#ifndef _PETPVC_HPP_
#define _PETPVC_HPP_

#include <iostream>
#include <memory>
#include <numeric>
#include <cassert>

#include "itkAddImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkThresholdImageFilter.h"

//TODO: Create unit tests for package.

namespace petpvc {

typedef itk::Image<float, 3> ImageType3D;
typedef itk::Image<float, 4> ImageType4D;

typedef itk::Image<short, 3> MaskImageType3D;

typedef itk::BinaryThresholdImageFilter<MaskImageType3D, MaskImageType3D> BinaryThresholdImageFilterType;
typedef itk::DiscreteGaussianImageFilter<ImageType3D, ImageType3D> BlurringImageFilterType;

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
              float val,
              typename TImageType::Pointer outputImage)
{
  typedef typename itk::MultiplyImageFilter<TImageType> FilterType;
  typename FilterType::Pointer multiply = FilterType::New();
  multiply->SetInput(a);
  multiply->SetConstant(val);
  multiply->Update();

  outputImage->Graft(multiply->GetOutput());
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
  typedef typename itk::DivideImageFilter<TImageType, TImageType, TImageType> FilterType;
  typename FilterType::Pointer divide = FilterType::New();
  divide->SetInput1(a);
  divide->SetInput2(b);
  divide->Update();

  outputImage->Graft(divide->GetOutput());
}

template<typename TInputImage>
void Duplicate(const typename TInputImage::Pointer input,
               typename TInputImage::Pointer &output){

  typedef itk::ImageDuplicator< TInputImage > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(input);
  duplicator->Update();
  output = duplicator->GetOutput();

}

template<typename TImageType=ImageType3D>
void RemoveNegativeNumbers(const typename TImageType::Pointer a,
                           typename TImageType::Pointer outputImage) {

  typedef typename itk::ThresholdImageFilter<TImageType> FilterType;
  typename FilterType::Pointer thresholdFilter = FilterType::New();

  thresholdFilter->ThresholdBelow( 0.0f );
  thresholdFilter->SetOutsideValue( 0.0f );
  thresholdFilter->SetInput(a);
  thresholdFilter->Update();

  outputImage->Graft(thresholdFilter->GetOutput());
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

  if (n < 0) {
    std::cerr << "Requested region " << n << " out of range!" << std::endl;
    throw false;
  }

  desiredStart[3] = n;
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

void GetRegion(const petpvc::ImageType4D::Pointer inImage, const int n, petpvc::ImageType3D::Pointer &outImage){
  petpvc::Extract4DVolTo3D<ImageType4D, ImageType3D>(inImage, n, outImage);
}

void GetRegion(petpvc::MaskImageType3D::Pointer inImage, const int n, petpvc::MaskImageType3D::Pointer &outImage) {

  typename BinaryThresholdImageFilterType::Pointer binThreshFilter = BinaryThresholdImageFilterType::New();

  const int targetIdx = n;

  //std::cout << "Extracting constant region " << targetIdx << std::endl;

  binThreshFilter->SetInput( inImage );
  binThreshFilter->SetInsideValue(1);
  binThreshFilter->SetOutsideValue(0);

  binThreshFilter->SetLowerThreshold( targetIdx );
  binThreshFilter->SetUpperThreshold( targetIdx );
  binThreshFilter->Update();

  outImage->Graft(binThreshFilter->GetOutput());

}

void GetRegionIndexList(const typename MaskImageType3D::Pointer input, std::vector<int> &indexList){

  //Get labels
  //For getting information about the mask labels
  typedef itk::LabelStatisticsImageFilter<MaskImageType3D, MaskImageType3D> LabelStatisticsFilterType;
  typedef typename LabelStatisticsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  typedef typename LabelStatisticsFilterType::LabelPixelType LabelPixelType;

  typename LabelStatisticsFilterType::Pointer labelStatsFilter = LabelStatisticsFilterType::New();
  labelStatsFilter->SetInput(input);
  labelStatsFilter->SetLabelInput(input);
  labelStatsFilter->Update();

  const int numOfRegions = labelStatsFilter->GetNumberOfLabels()-1;

  std::cout << "Number of labels found: " << numOfRegions << std::endl;

  if ( numOfRegions == 0) {
    std::cerr << "No regions found!" << std::endl;
    return;
  }

  //Put labels into array/vector
  std::vector<int> v;
  v.reserve(numOfRegions+1);

  std::cout << "Found Label(s): ";
  for( typename ValidLabelValuesType::const_iterator vIt=labelStatsFilter->GetValidLabelValues().begin();
       vIt != labelStatsFilter->GetValidLabelValues().end(); ++vIt) {
    if (labelStatsFilter->HasLabel(*vIt)) {
      LabelPixelType labelValue = *vIt;
      std::cout << labelValue << " ";
      v.push_back(labelValue);
    }
  }

  std::cout << std::endl;
  std::cout << "Label indices: " ;

  for (auto x : v) {
    std::cout <<  x << " ";
  }
  std::cout << std::endl;
  indexList = v;

}

void GetRegionIndexList(const typename ImageType4D::Pointer input, std::vector<int> &indexList) {

  const typename ImageType4D::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();

  std::cout << "Dimension = " << inputSize.Dimension << std::endl;

  if (inputSize.Dimension == 4){

    std::vector<int> v(inputSize[3]);
    std::iota(std::begin(v), std::end(v), 0);
    indexList = v;
    return;
  }
}

void GetRegionalMeans(const ImageType3D::Pointer input, const MaskImageType3D::Pointer mask,
                      std::vector<float> &meansList){

  //For getting information about the mask labels
  typedef itk::LabelStatisticsImageFilter<ImageType3D, MaskImageType3D> LabelStatisticsFilterType;
  typedef typename LabelStatisticsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  typedef typename LabelStatisticsFilterType::LabelPixelType LabelPixelType;

  typename LabelStatisticsFilterType::Pointer labelStatsFilter = LabelStatisticsFilterType::New();
  labelStatsFilter->SetInput(input);
  labelStatsFilter->SetLabelInput(mask);
  labelStatsFilter->Update();

  const int numOfRegions = labelStatsFilter->GetNumberOfLabels()-1;

  //Put means into array/vector
  std::vector<float> v;
  v.reserve(numOfRegions+1);

  std::cout << "Found 3D means(s): ";
  for( typename ValidLabelValuesType::const_iterator vIt=labelStatsFilter->GetValidLabelValues().begin();
       vIt != labelStatsFilter->GetValidLabelValues().end(); ++vIt) {
    if (labelStatsFilter->HasLabel(*vIt)) {
      LabelPixelType labelValue = *vIt;
      float fNewRegMean = std::max( labelStatsFilter->GetMean( labelValue ), 0.0 );
      std::cout << fNewRegMean << " ";
      v.push_back(fNewRegMean);
    }
  }

  std::cout << std::endl;
  meansList = v;
  /*
  std::cout << "meansList = ";
  for (auto x : meansList){
    std::cout << x << " ";
  }
  std::cout << std::endl;*/

}

void GetRegionalMeans(const ImageType3D::Pointer input, const ImageType4D::Pointer mask,
                        std::vector<float> &meansList){

  std::cout << "Calculating 4D regional means" << std::endl;

  const typename ImageType4D::SizeType maskSize = mask->GetLargestPossibleRegion().GetSize();
  std::cout << "Dimension of mask = " << maskSize.Dimension << std::endl;

  int numOfRegions = maskSize[3];

  std::vector<float> v(numOfRegions);

  ImageType3D::Pointer tmpMask = ImageType3D::New();
  ImageType3D::Pointer clippedPET = ImageType3D::New();

  typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
  StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();

  for (int n=0; n < numOfRegions; n++){
    ImageType3D::Pointer tmpMask = ImageType3D::New();
    Extract4DVolTo3D<ImageType4D, ImageType3D>(mask,n,tmpMask);

    statsFilter->SetInput(tmpMask);
    statsFilter->Update();

    //Get sum of the clipped image.
    float sumOfMaskReg = statsFilter->GetSum();

    Multiply(input,tmpMask,clippedPET);

    statsFilter->SetInput(clippedPET);
    statsFilter->Update();

    //Get sum of the clipped image.
    float sumOfPETReg = statsFilter->GetSum();

    //Place regional mean into vector.
    float regMean = std::max( sumOfPETReg / sumOfMaskReg, 0.0f);
    v[n] = regMean;
    //std::cout << std::endl << "Sum = " << fSumOfPETReg << " , " << "Mean = " << vecRegMeansCurrent.get( i-1 ) << " , Size = " << vecRegSize.get(i - 1) << std::endl;

  }
  meansList = v;
  std::cout << "Total no. of means: " << meansList.size() << std::endl;

  std::cout << "Means : ";

  for ( auto i : meansList){
    std::cout << i << " ";
  }

  std::cout << std::endl;
}

void GetVolume(const petpvc::ImageType4D::Pointer inImage, const int n, petpvc::ImageType3D::Pointer &outImage){
  petpvc::Extract4DVolTo3D<ImageType4D, ImageType3D>(inImage, n, outImage);
}

void GetVolume(const petpvc::ImageType3D::Pointer inImage, const int n, petpvc::ImageType3D::Pointer &outImage){

  if (n != 0) {
    std::cerr << "Volume does not exist!" << std::endl;
    std::cout << "Returning 1st volume" << std::endl;
  }

  Duplicate<ImageType3D>(inImage,outImage);
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

//TODO: Implement paste into 4D
void PasteInto(const ImageType3D::Pointer inVol,
               const int n,
               ImageType4D::Pointer destVol) {

  typedef itk::CastImageFilter< ImageType3D, ImageType4D > CastImageFilterType;
  typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
  castFilter->SetInput( inVol );
  castFilter->Update();

  typedef typename itk::PasteImageFilter<ImageType4D> PasteFilterType;
  typename PasteFilterType::Pointer pasteFilter = PasteFilterType::New();

  typename ImageType4D::IndexType desiredStart;
  desiredStart.Fill(0);

  typename ImageType4D::SizeType desiredSize =
      destVol->GetLargestPossibleRegion().GetSize();

  if (n < 0) {
    std::cerr << "Requested region " << n << " out of range!" << std::endl;
    throw false;
  }

  desiredStart[3] = n;

  pasteFilter->SetSourceImage( castFilter->GetOutput() );
  pasteFilter->SetDestinationImage( destVol );
  pasteFilter->SetDestinationIndex( desiredStart );
  pasteFilter->SetSourceRegion( castFilter->GetOutput()->GetLargestPossibleRegion() );
  pasteFilter->Update();

  destVol->Graft(pasteFilter->GetOutput());

}

template<typename T>
int GetNumberOfVolumes(const typename T::Pointer img ){
  const int dims = T::ImageDimension;

  if (dims < 3){
    std::cerr << "Must be 3D or 4D image" << std::endl;
    return -1;
  }

  if (dims == 3)
    return 1;

  if (dims == 4)
    return img->GetLargestPossibleRegion().GetSize(3);

  return -1;
}

void ConfigureGaussian(typename BlurringImageFilterType::Pointer &gaussian,
                        const float x_in_mm, float y_in_mm = -1.0f, float z_in_mm = -1.0f){

  if (y_in_mm == -1.0f)
    y_in_mm = x_in_mm;

  if (z_in_mm == -1.0f)
    z_in_mm = x_in_mm;

  //Calculate the variance for a given FWHM.
  double psfVar[3];

  psfVar[0] = pow( x_in_mm/(2.0 * sqrt(2.0 * log(2.0))), 2.0);
  psfVar[1] = pow( y_in_mm/(2.0 * sqrt(2.0 * log(2.0))), 2.0);
  psfVar[2] = pow( z_in_mm/(2.0 * sqrt(2.0 * log(2.0))), 2.0);

  gaussian->SetVariance(psfVar);

}

template<typename TInputImage, typename TBlurFilter>
void ApplySmoothing(const typename TInputImage::Pointer pet,
                      const typename TBlurFilter::Pointer blur,
                      typename ImageType3D::Pointer &output){

  typedef itk::CastImageFilter< TInputImage, ImageType3D > CastImageFilterType;
  typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
  castFilter->SetInput( pet );
  castFilter->Update();

  blur->SetInput( castFilter->GetOutput() );
  blur->Update();

  Duplicate<ImageType3D>(blur->GetOutput(),output);

}

template<typename TInputImage=ImageType3D>
float CalculateSum(const typename TInputImage::Pointer input) {

  itk::ImageRegionIterator<TInputImage> it(input, input->GetLargestPossibleRegion());

  float imageSum = 0.0f;

  it.GoToBegin();

  while (!it.IsAtEnd()){
    imageSum += it.Get();
    ++it;
  }

  return imageSum;

}

template<typename TInputImage=ImageType3D>
float CalculateSumOfSquares(const typename TInputImage::Pointer a) {

  itk::ImageRegionIterator<TInputImage> it1(a, a->GetLargestPossibleRegion());
  float imageSumSq = 0.0f;

  it1.GoToBegin();

  float val;

  while (!it1.IsAtEnd()){
    val = it1.Get();
    imageSumSq += val*val;
    ++it1;
  }

  return imageSumSq;
}

template<typename TInputImage=ImageType3D>
float CalculateSumOfSquareDiff(const typename TInputImage::Pointer a,
                               const typename TInputImage::Pointer b) {

  itk::ImageRegionIterator<TInputImage> it1(a, a->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TInputImage> it2(b, b->GetLargestPossibleRegion());

  float imageSumSqDiff = 0.0f;

  it1.GoToBegin(); it2.GoToBegin();

  float diff;

  while (!it1.IsAtEnd()){
    diff = it1.Get() - it2.Get();
    imageSumSqDiff += diff*diff;
    ++it1; ++it2;
  }

  return imageSumSqDiff;

}

}; //end namespace petpvc

#endif