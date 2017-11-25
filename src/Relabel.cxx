/*
   Relabel.cxx

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

   Converts one 3-D parcellation into another.

*/

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkImageDuplicator.h>
#include <itkOneWayEquivalencyTable.h>
#include <metaCommand.h>

#include "EnvironmentInfo.h"

#include <iostream>
#include <fstream>
#include <map> 

typedef itk::Image<short, 3>   ImageType;

typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;

typedef itk::ImageRegionIterator<ImageType> ImageIterator;
typedef itk::ImageRegionConstIterator<ImageType> ConstImageIterator;

int main(int argc, char *argv[])
{

	const char * const AUTHOR = "Benjamin A. Thomas";
  const char * const APP_TITLE = "Relabel an image";

  std::stringstream version_number;
  version_number << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
  const std::string VERSION_NO = version_number.str();
	MetaCommand command;

	command.SetVersion( VERSION_NO.c_str() );
	command.SetAuthor( AUTHOR );
	command.SetName( APP_TITLE );
	command.SetDescription("Relabels a parcellation image");
	command.SetAcknowledgments( " " );
	command.SetCategory("PETPVC");

	command.SetOption("Input", "i", true, "Input file");
  command.AddOptionField("Input", "infilename", MetaCommand::FILE, true, "", "", MetaCommand::DATA_IN);

  command.SetOption("Output", "o", true, "Output file");
  command.AddOptionField("Output", "outfilename", MetaCommand::FILE, true, "", "", MetaCommand::DATA_OUT);

  command.SetOption("Parcellation", "p", true, "Description file");
	command.SetOptionLongTag("Parcellation", "parc");
	command.AddOptionField("Parcellation", "parcfile", MetaCommand::FILE, true, "", "", MetaCommand::DATA_IN);

  command.SetOption("Type", "t", true, "Parcellation type");
	command.SetOptionLongTag("Type", "type");
  command.AddOptionField("Type", "parctype", MetaCommand::STRING, true, "", "");

	if( !command.Parse(argc,argv) )
	{
		return EXIT_FAILURE;
	}	

	std::string inFileName = command.GetValueAsString("Input", "infilename");
	std::string maskDescriptionFileName = command.GetValueAsString("Parcellation", "parcfile");
	std::string outputFileName = command.GetValueAsString("Output", "outfilename");
	std::string targetColumnName = command.GetValueAsString("Type", "parctype");

	std::cout << "Input file: " << inFileName << std::endl;
	std::cout << "Output file: " << outputFileName << std::endl;
	std::cout << "Description file: " << maskDescriptionFileName << std::endl;
	std::cout << "Parcellation type: " << targetColumnName << std::endl;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inFileName );

	std::ifstream maskDesriptionFile;
	maskDesriptionFile.open( maskDescriptionFileName.c_str() , std::ios::in );

	std::vector<std::string> columnHeaders;
	std::string line;

	std::map<float, float> eqTable;

	if ( maskDesriptionFile.is_open() ){

		int targetColumnIndex =0;

		//Get header of CSV
		getline( maskDesriptionFile, line );

		//Remove \r
		if ( line.size() && line[line.size()-1] == '\r' ){
           line = line.substr( 0, line.size() - 1 );
		}

		std::istringstream s(line);
		std::string field;

		//Find requested parcellation type.
		int count = 0;
		while (getline(s, field,',')) {
			//std::cout << count << " : " << field << std::endl;
			columnHeaders.push_back(field);

			if (field == targetColumnName){
				targetColumnIndex = count;
			}

			count++;
		}

		if ( targetColumnIndex != 0) {
			std::cout << "Found desired parcellation index in column : " << targetColumnIndex << std::endl << std::endl;
		}
		else {
			std::cerr << "[Error]\tCannot find desired parcellation " << targetColumnName << 
				" in " << maskDescriptionFileName << "!" << std::endl;

			return EXIT_FAILURE;

		}

		std::string currentRegion;
		float sourceVal;
		float destVal;

		std::cout << "Mapping:" << std::endl;

		//Build equivalency table mapping.
		while ( getline( maskDesriptionFile, line ) ) {

			std::istringstream ss(line);
			std::vector<std::string> rowValues;
			
			while (getline(ss, field,',')) {
				rowValues.push_back(field);
			}	

			currentRegion = rowValues[0];

			std::stringstream inval;
			inval << rowValues[1];
			inval >> sourceVal;

			inval.clear();
			inval.str("");
			inval << rowValues[targetColumnIndex];
			inval >> destVal;

			if ( destVal != 0) {
				eqTable[sourceVal]=destVal;
				std::cout << "\tID: "	<< sourceVal << " -> " << eqTable[sourceVal] 
					<<  "\t\t(" << currentRegion << ")" << std::endl;
			}
		}
	}
	else {
		std::cerr << "[Error]\tCannot open mask description file: " << maskDescriptionFileName << "!" << std::endl;
		return EXIT_FAILURE;
	}

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cerr << "[Error]\tImage file: " << inFileName << " cannot be loaded!" << std::endl;
		return EXIT_FAILURE;
	}

	typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  	DuplicatorType::Pointer duplicator = DuplicatorType::New();
  	duplicator->SetInputImage( reader->GetOutput() );
  	duplicator->Update();

  	ImageType::Pointer outputImage = duplicator->GetOutput();

	ImageIterator outputIt( outputImage, outputImage->GetRequestedRegion() );
	
	while (! outputIt.IsAtEnd() )
	{
		outputIt.Set( 0 );
		++outputIt;
	}	

	outputIt.GoToBegin();

	ConstImageIterator inputIt( reader->GetOutput(), reader->GetOutput()->GetRequestedRegion() );

	//Apply eqivalency table 
	while ( !inputIt.IsAtEnd() )
	{
		if ( eqTable.find( inputIt.Get() ) != eqTable.end() ) {
			outputIt.Set( eqTable[inputIt.Get()] );
		}
		++inputIt;
		++outputIt;
	}	

	//Write to disk.
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outputFileName );
	writer->SetInput( outputImage );
	
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
