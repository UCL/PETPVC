/*
   Make4D.cxx

   Author:      Benjamin A. Thomas

   Copyright 2017 Institute of Nuclear Medicine, University College London.
   Copyright 2014 Clinical Imaging Research Centre, A*STAR-NUS.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Takes a 3-D mask image and creates a 4-D image with 1 volume per label.

*/

#include <iostream>
#include <fstream>
#include <map>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageDuplicator.h>
#include <itkLabelGeometryImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkJoinSeriesImageFilter.h>
#include <metaCommand.h>

#include "EnvironmentInfo.h"

typedef itk::Image<short, 3> ImageType;
typedef itk::Image<short, 3> MaskImageType;
typedef itk::Image<short, 4> ImageType4D;

typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileWriter<ImageType4D> WriterType;

typedef itk::ImageRegionIterator<ImageType> ImageIterator;
typedef itk::ImageRegionConstIterator<ImageType> ConstImageIterator;

typedef itk::LabelGeometryImageFilter<ImageType> LabelGeometryImageFilterType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryThresholdImageFilterType;

typedef itk::JoinSeriesImageFilter<ImageType, ImageType4D> JoinSeriesImageFilterType;

int main(int argc, char *argv[])
{

    const char * const AUTHOR = "Benjamin A. Thomas";
    const char * const APP_TITLE = "3-D to 4-D mask creation";

    std::stringstream version_number;
    version_number << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
    const std::string VERSION_NO = version_number.str();

	MetaCommand command;

	command.DisableDeprecatedWarnings();

	command.SetVersion(VERSION_NO.c_str());
	command.SetAuthor(AUTHOR);
	command.SetName(APP_TITLE);
	command.SetDescription("Creates a 4-D region mask image from a 3-D mask image");
	command.SetAcknowledgments(" ");
	command.SetCategory("PETPVC");

	command.SetOption("Input", "i", true, "Input file");
	command.AddOptionField("Input", "infilename", MetaCommand::FILE, true, "", "", MetaCommand::DATA_IN);

	command.SetOption("Output", "o", true, "Output file");
	command.AddOptionField("Output", "outfilename", MetaCommand::FILE, true, "", "", MetaCommand::DATA_OUT);

	if (!command.Parse(argc, argv))
	{
		return EXIT_FAILURE;
	}

	std::string inFileName = command.GetValueAsString("Input", "infilename");
	std::string outputFileName = command.GetValueAsString("Output", "outfilename");

	std::cout << "Input file: " << inFileName << std::endl;
	std::cout << "Output file: " << outputFileName << std::endl;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inFileName);

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cerr << "[Error]\tImage file: " << inFileName << " cannot be loaded!" << std::endl;
		return EXIT_FAILURE;
	}

	ImageType::Pointer inputImage = reader->GetOutput();

	LabelGeometryImageFilterType::Pointer labelGeometryFilter = LabelGeometryImageFilterType::New();
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput(inputImage);
	labelGeometryImageFilter->SetIntensityInput(inputImage);

	labelGeometryImageFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels =
	labelGeometryImageFilter->GetLabels();

	std::cout << "No. of labels found:\t" << allLabels.size() << std::endl;

	BinaryThresholdImageFilterType::Pointer binThreshFilter = BinaryThresholdImageFilterType::New();
	binThreshFilter->SetInput(inputImage);
	binThreshFilter->SetInsideValue(1);
	binThreshFilter->SetOutsideValue(0);

	JoinSeriesImageFilterType::Pointer joinFilter = JoinSeriesImageFilterType::New();
	ImageType::Pointer maskImage = ImageType::New();

	std::cout << "Mapping:" << std::endl;

	for (int i = 0; i < allLabels.size(); i++)
	{
		std::cout << "\tID: " << allLabels[i] << " -> volume " << i+1 << std::endl;
		int labelValue = allLabels[i];
		binThreshFilter->SetLowerThreshold(labelValue);
		binThreshFilter->SetUpperThreshold(labelValue);
		binThreshFilter->Update();

		maskImage = binThreshFilter->GetOutput();

		joinFilter->SetInput(i, maskImage);
		maskImage->DisconnectPipeline();
	}

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputFileName);
	writer->SetInput(joinFilter->GetOutput());

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cerr << "[Error]\tCannot write output image: " << outputFileName << "!" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
