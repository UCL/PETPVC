/*
   VanCittert.cxx

   Author:      Benjamin A. Thomas

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

   This program implements the reblurred Van-Cittert (RVC) partial volume correction
   technique. Please cite the following paper:

	Tohka, J. and Reilhac A., (2008). "Deconvolution-based partial volume
	correction in Raclopride-PET and Monte Carlo comparison to MR-based
	method", NeuroImage, vol. 39. 1570--1584.        

 */

#include "EnvironmentInfo.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <metaCommand.h>

#include "petpvcVanCittertPVCImageFilter.h"

typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<float, 3> PETImageType;

typedef itk::ImageFileReader<PETImageType> PETReaderType;
typedef itk::ImageFileWriter<PETImageType> PETWriterType;

//Produces the text for the acknowledgments dialog in Slicer.
std::string getAcknowledgments(void);

int main(int argc, char *argv[])
{
    const char * const AUTHOR = "Benjamin A. Thomas";
    const char * const APP_TITLE = "Reblurred Van-Cittert (RVC)";

    std::stringstream version_number;
    version_number << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
    const std::string VERSION_NO = version_number.str();

    typedef petpvc::VanCittertPVCImageFilter< PETImageType >  FilterType;

    //Setting up command line argument list.
    MetaCommand command;

    command.SetVersion(VERSION_NO.c_str());
    command.SetAuthor(AUTHOR);
    command.SetName(APP_TITLE);
    command.SetDescription(
        "Performs reblurred Van-Cittert (RVC) partial volume correction");

    std::string sAcks = getAcknowledgments();
    command.SetAcknowledgments(sAcks.c_str());

    command.SetCategory("PETPVC");

    command.AddField("petfile", "PET filename", MetaCommand::IMAGE, MetaCommand::DATA_IN);
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

    command.SetOption("Iterations", "i", false, "Number of iterations");
    command.SetOptionLongTag("Iterations", "iter");
    command.AddOptionField("Iterations", "Val", MetaCommand::INT, false, "30");

    command.SetOption("Alpha", "a", false, "Alpha value");
    command.SetOptionLongTag("Alpha", "alpha");
    command.AddOptionField("Alpha", "aval", MetaCommand::FLOAT, false, "1.5");

    command.SetOption("Stop", "s", false, "Stopping criterion");
    command.SetOptionLongTag("Stop", "stop");
    command.AddOptionField("Stop", "stopval", MetaCommand::FLOAT, false, "0.01");

    command.SetOption("debug", "d", false,"Prints debug information");
    command.SetOptionLongTag("debug", "debug");

    //Parse command line.
    if (!command.Parse(argc, argv)) {
        return EXIT_FAILURE;
    }

    //Get image filenames
    std::string sPETFileName = command.GetValueAsString("petfile");
    std::string sOutputFileName = command.GetValueAsString("outputfile");

    //Get values for PSF.
    float fFWHM_x = command.GetValueAsFloat("FWHMx", "X");
    float fFWHM_y = command.GetValueAsFloat("FWHMy", "Y");
    float fFWHM_z = command.GetValueAsFloat("FWHMz", "Z");

    //Get number of iterations
    int nNumOfIters = command.GetValueAsInt("Iterations", "Val");

    //Get value for alpha.
    float fAlpha = command.GetValueAsFloat("Alpha", "aval");

    //Get value for stopping criterion.
    float fStop = command.GetValueAsFloat("Stop", "stopval");
	
    //Make vector of FWHM in x,y and z.
    VectorType vFWHM;
    vFWHM[0] = fFWHM_x;
    vFWHM[1] = fFWHM_y;
    vFWHM[2] = fFWHM_z;

    //Toggle debug mode
    bool bDebug = command.GetValueAsBool("debug");

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

    FilterType::Pointer vcFilter = FilterType::New();
    vcFilter->SetInput( petReader->GetOutput() );
    vcFilter->SetPSF(vVariance);
    vcFilter->SetIterations( nNumOfIters );
	vcFilter->SetAlpha( fAlpha );
    vcFilter->SetStoppingCond( fStop );
    vcFilter->SetVerbose ( bDebug );

    //Perform VC.
    try {
        vcFilter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tfailure applying Van Cittert on: " << sPETFileName
                  << "\n" << err
                  << std::endl;
        return EXIT_FAILURE;
    }

    PETWriterType::Pointer petWriter = PETWriterType::New();
    petWriter->SetFileName(sOutputFileName);
    petWriter->SetInput( vcFilter->GetOutput() );

    try {
        petWriter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "\n[Error]\tCannot write output file: " << sOutputFileName
                  << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

std::string getAcknowledgments(void)
{
    //Produces acknowledgments string for 3DSlicer.
    std::string sAck = "This program implements the reblurred Van-Cittert (RVC) partial volume correction technique. Please cite the following paper:\n"
                       "\tTohka, J. and Reilhac A., (2008). \"Deconvolution-based partial volume correction in Raclopride-PET\n\tand Monte Carlo comparison to MR-based method\", NeuroImage, vol. 39. 1570--1584.";

    return sAck;
}

