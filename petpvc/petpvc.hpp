#pragma once

#ifndef _PETPVC_HPP
#define _PETPVC_HPP

#include <iostream>
#include <memory>
#include <numeric>

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
#include "itkSubtractImageFilter.h"

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

void GetRegion(const petpvc::ImageType4D::Pointer inImage, const int n, petpvc::ImageType3D::Pointer &outImage){
  petpvc::Extract4DVolTo3D<ImageType4D, ImageType3D>(inImage, n, outImage);
}

void GetRegion(petpvc::MaskImageType3D::Pointer inImage, const int n, petpvc::MaskImageType3D::Pointer &outImage) {

  typename BinaryThresholdImageFilterType::Pointer binThreshFilter = BinaryThresholdImageFilterType::New();

  const int targetIdx = n;

  std::cout << "Extracting region " << targetIdx << std::endl;

  binThreshFilter->SetInput( inImage );
  binThreshFilter->SetInsideValue(1);
  binThreshFilter->SetOutsideValue(0);

  binThreshFilter->SetLowerThreshold( targetIdx );
  binThreshFilter->SetUpperThreshold( targetIdx );
  binThreshFilter->Update();

  outImage->Graft(binThreshFilter->GetOutput());

}

template<typename TImageType>
void GetRegionIndexList(const typename TImageType::Pointer input, std::vector<int> &indexList){

  const typename TImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();

  if (inputSize.Dimension == 4){

    std::vector<int> v(inputSize[3]);
    std::iota(std::begin(v), std::end(v), 0);
    indexList = v;
    return;
  }

  if (inputSize.Dimension == 3){
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

}

void GetVolume(const petpvc::ImageType4D::Pointer inImage, const int n, petpvc::ImageType3D::Pointer &outImage){
  petpvc::Extract4DVolTo3D<ImageType4D, ImageType3D>(inImage, n, outImage);
}

void GetVolume(const petpvc::ImageType3D::Pointer inImage, const int n, petpvc::ImageType3D::Pointer &outImage){

  if (n != 0) {
    std::cerr << "Volume does not exist!" << std::endl;
    std::cout << "Returning 1st volume" << std::endl;
  }

  outImage->Graft( inImage );
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

  if (n - 1 < 0) {
    std::cerr << "Requested region " << n << " out of range!" << std::endl;
    throw false;
  }

  desiredStart[3] = n - 1;
  desiredSize[3] = 0;

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
                        const float x_in_mm, float y_in_mm=0, float z_in_mm=0){

  if (y_in_mm == 0)
    y_in_mm = x_in_mm;

  if (z_in_mm == 0)
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
                      typename TInputImage::Pointer &output){

    blur->SetInput( pet );
    blur->Update();

    output->Graft(blur->GetOutput());

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

}; //end namespace petpvc

#endif