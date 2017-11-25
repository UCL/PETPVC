/*
   MullerGartner.cxx

   Author:      Benjamin A. Thomas

   Copyright 2013-2015 Institute of Nuclear Medicine, University College London.
   Copyright 2015 Clinical Imaging Research Centre, A*STAR-NUS.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   This program implements the Muller-Gartner (MG) partial volume
   correction (PVC) technique. The method is described in:
        Muller-Gartner, H. W. et al. (1992). "Measurement of radiotracer
        concentration in brain gray matter using positron emission
        tomography: MRI-based correction for partial volume effects".
        J Cereb Blood Flow Metab, 12(4), 571-83.

 */

#include "EnvironmentInfo.h"
#include <string>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include "petpvcMullerGartnerImageFilter.h"

#include <metaCommand.h>

typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<float, 4> MaskImageType;
typedef itk::Image<float, 3> PETImageType;

typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

typedef itk::ImageFileReader<PETImageType> PETReaderType;
typedef itk::ImageFileWriter<PETImageType> PETWriterType;

//Extracts a 3D volume from 4D file.
typedef itk::ExtractImageFilter<MaskImageType, PETImageType> ExtractFilterType;

typedef itk::ImageDuplicator<PETImageType> DuplicatorType;

typedef petpvc::MullerGartnerImageFilter< PETImageType, PETImageType, PETImageType, PETImageType > MGFilterType;
//Function definitions:

//Produces the text for the acknowledgments dialog in Slicer.
std::string getAcknowledgments(void);

int main(int argc, char *argv[])
{
    const char * const AUTHOR = "Benjamin A. Thomas";
    const char * const APP_TITLE = "Muller-Gartner (MG) PVC";

    std::stringstream version_number;
    version_number << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
    const std::string VERSION_NO = version_number.str();

    //Setting up command line argument list.
    MetaCommand command;

    command.SetVersion(VERSION_NO.c_str());
    command.SetAuthor(AUTHOR);
    command.SetName(APP_TITLE);
    command.SetDescription(
        "Performs Muller-Gartner (MG) partial volume correction");

    std::string sAcks = getAcknowledgments();
    command.SetAcknowledgments(sAcks.c_str());

    command.SetCategory("PETPVC");

    command.AddField("petfile", "PET filename", MetaCommand::IMAGE, MetaCommand::DATA_IN);
    command.AddField("maskfile", "mask filename", MetaCommand::IMAGE, MetaCommand::DATA_IN);
    command.AddField("outputfile", "output filename", MetaCommand::IMAGE, MetaCommand::DATA_OUT);

    command.SetOption("FWHMx", "x", true,
                      "The full-width at half maximum in mm along x-axis");
    command.AddOptionField("FWHMx", "X", MetaCommand::FLOAT, true, "");

    command.SetOption("FWHMy", "y", true,
                      "The full-width at half maximum in mm along y-axis");
    command.AddOptionField("FWHMy", "Y", MetaCommand::FLOAT, true, "");

    command.SetOption("FWHMz", "z", true,
                      "The full-width at half maximum in mm along z-axis");
    command.AddOptionField("FWHMz", "Z", MetaCommand::FLOAT, true, "");

    command.SetOption("debug", "d", false,"Prints debug information");
    command.SetOptionLongTag("debug", "debug");

    //Parse command line.
    if (!command.Parse(argc, argv)) {
        return EXIT_FAILURE;
    }

    //Get image filenames
    std::string sPETFileName = command.GetValueAsString("petfile");
    std::string sMaskFileName = command.GetValueAsString("maskfile");
    std::string sOutputFileName = command.GetValueAsString("outputfile");

    //Get values for PSF.
    float fFWHM_x = command.GetValueAsFloat("FWHMx", "X");
    float fFWHM_y = command.GetValueAsFloat("FWHMy", "Y");
    float fFWHM_z = command.GetValueAsFloat("FWHMz", "Z");

    //Make vector of FWHM in x,y and z.
    VectorType vFWHM;
    vFWHM[0] = fFWHM_x;
    vFWHM[1] = fFWHM_y;
    vFWHM[2] = fFWHM_z;

    //Toggle debug mode
    bool bDebug = command.GetValueAsBool("debug");

    //Create reader for mask image.
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName(sMaskFileName);

    //Try to read mask.
    try {
        maskReader->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tCannot read mask input file: " << sMaskFileName
                  << std::endl;
        return EXIT_FAILURE;
    }

    //Create reader for PET image.
    PETReaderType::Pointer petReader = PETReaderType::New();
    petReader->SetFileName(sPETFileName);

    //Try to read PET.
    try {
        petReader->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tCannot read PET input file: " << sPETFileName
                  << std::endl;
        return EXIT_FAILURE;
    }

    //Calculate the variance for a given FWHM.
    VectorType vVariance;
    vVariance = vFWHM / (2.0 * sqrt(2.0 * log(2.0)));
    //std::cout << vVariance << std::endl;

    VectorType vVoxelSize = petReader->GetOutput()->GetSpacing();
    //std::cout << vVoxelSize << std::endl;

    vVariance[0] = pow(vVariance[0], 2);
    vVariance[1] = pow(vVariance[1], 2);
    vVariance[2] = pow(vVariance[2], 2);

    //Get mask image size.
    MaskImageType::SizeType imageSize =
        maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();

    //If mask is not 4D, then quit.
    if (imageSize.Dimension != 4) {
        std::cerr << "[Error]\tMask file: " << sMaskFileName << " must be 4-D!"
                  << std::endl;
        return EXIT_FAILURE;
    }

    MaskImageType::IndexType desiredStart;
    desiredStart.Fill(0);
    MaskImageType::SizeType desiredSize = imageSize;

    //Extract filter used to extract 3D volume from 4D file.
    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput(maskReader->GetOutput());
    extractFilter->SetDirectionCollapseToIdentity(); // This is required.

    PETImageType::Pointer imageGM = PETImageType::New();
    PETImageType::Pointer imageWM = PETImageType::New();
    //PETImageType::Pointer imageCSF = PETImageType::New();

    //Starts reading from 4D volume at index (0,0,0,i) through to
    //(maxX, maxY, maxZ,0), i.e. one 3D brain mask.
    desiredStart[3] = 0;
    desiredSize[3] = 0;

    //Get GM mask.
    extractFilter->SetExtractionRegion(
        MaskImageType::RegionType(desiredStart, desiredSize));
    extractFilter->Update();

    imageGM = extractFilter->GetOutput();
    imageGM->SetDirection(petReader->GetOutput()->GetDirection());
    imageGM->UpdateOutputData();
    imageGM->DisconnectPipeline();

    desiredStart[3] = 1;

    //Get WM mask.
    extractFilter->SetExtractionRegion(
        MaskImageType::RegionType(desiredStart, desiredSize));
    extractFilter->Update();

    imageWM = extractFilter->GetOutput();
    imageWM->SetDirection(petReader->GetOutput()->GetDirection());
    imageWM->UpdateOutputData();
    imageWM->DisconnectPipeline();

    desiredStart[3] = 2;

    //Get CSF mask.
    /*
    extractFilter->SetExtractionRegion(
        MaskImageType::RegionType(desiredStart, desiredSize));
    extractFilter->Update();

    imageCSF = extractFilter->GetOutput();
    imageCSF->SetDirection(petReader->GetOutput()->GetDirection());
    imageCSF->UpdateOutputData();
    imageCSF->DisconnectPipeline();*/

    //Set-up Muller-Gartner filter
    MGFilterType::Pointer MGFilter = MGFilterType::New();

    MGFilter->SetInput1(petReader->GetOutput());
    MGFilter->SetInput2(imageGM);
    MGFilter->SetInput3(imageWM);
    MGFilter->SetVerbose( bDebug );
    MGFilter->SetPSF(vVariance);
    MGFilter->SetWM(0);

    //Set-up output file writer
    PETWriterType::Pointer petWriter = PETWriterType::New();

    petWriter->SetInput(MGFilter->GetOutput());
    petWriter->SetFileName(sOutputFileName.c_str());

    //Write file
    try {
        petWriter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tCannot write output file: " << sOutputFileName
                  << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

std::string getAcknowledgments(void)
{
    //Produces acknowledgments string for 3DSlicer.
    std::string sAck = "This program implements the Muller-Gartner (MG) partial volume correction (PVC) technique.\n"
                       "The method is described in:\n"
                       "\tMuller-Gartner, H. W. et al. (1992). \"Measurement of radiotracer\n"
                       "\tconcentration in brain gray matter using positron emission\n"
                       "\ttomography: MRI-based correction for partial volume effects.\"\n"
                       "\tJ Cereb Blood Flow Metab, 12(4), 571-83.";
    return sAck;
}
