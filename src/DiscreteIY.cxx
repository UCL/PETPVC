/*
   DiscreteIY.cxx

   Author:      Benjamin A. Thomas

   Copyright 2014 Institute of Nuclear Medicine, University College London.
   Copyright 2014 Clinical Imaging Research Centre, A*STAR-NUS

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   This program implements the Iterative Yang (IY) partial volume correction
   (PVC) technique. Please cite the following paper:
        Erlandsson, K. and Buvat, I. and Pretorius, P.H. and Thomas, B.A.
        and Hutton, B.F., (2012). "A review of partial volume correction
        techniques for emission tomography and their applications in neurology,
        cardiology and oncology", Physics in Medicine and Biology,
        vol. 57, no. 21, R119-59.

 */

#include "EnvironmentInfo.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <metaCommand.h>
#include <itkTimeProbe.h>

#include "petpvcDiscreteIYPVCImageFilter.h"

typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<short, 3> MaskImageType;
typedef itk::Image<float, 3> PETImageType;

typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
typedef itk::ImageFileReader<PETImageType> PETReaderType;
typedef itk::ImageFileWriter<PETImageType> PETWriterType;

//Produces the text for the acknowledgments dialog in Slicer.
std::string getAcknowledgments(void);

int main(int argc, char *argv[])
{
    const char * const AUTHOR = "Benjamin A. Thomas";
    const char * const APP_TITLE = "Discrete Iterative Yang (IY) PVC";

    std::stringstream version_number;
    version_number << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
    const std::string VERSION_NO = version_number.str();

    typedef petpvc::DiscreteIYPVCImageFilter<PETImageType, MaskImageType>  FilterType;

    //Setting up command line argument list.
    MetaCommand command;

    command.SetVersion(VERSION_NO.c_str());
    command.SetAuthor(AUTHOR);
    command.SetName(APP_TITLE);
    command.SetDescription(
        "Performs discrete iterative Yang (IY) partial volume correction");

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

    command.SetOption("Iterations", "i", false, "Number of iterations");
    command.SetOptionLongTag("Iterations", "iter");
    command.AddOptionField("Iterations", "Val", MetaCommand::INT, false, "10");

    command.SetOption("debug", "d", false,"Prints debug information");
    command.SetOptionLongTag("debug", "debug");

    //Parse command line.
    if (!command.Parse(argc, argv)) {
        return EXIT_FAILURE;
    }

    itk::TimeProbe timer;

    timer.Start();

    //Get image filenames
    std::string sPETFileName = command.GetValueAsString("petfile");
    std::string sMaskFileName = command.GetValueAsString("maskfile");
    std::string sOutputFileName = command.GetValueAsString("outputfile");

    //Get values for PSF.
    float fFWHM_x = command.GetValueAsFloat("FWHMx", "X");
    float fFWHM_y = command.GetValueAsFloat("FWHMy", "Y");
    float fFWHM_z = command.GetValueAsFloat("FWHMz", "Z");

    //Get number of iterations
    int nNumOfIters = command.GetValueAsInt("Iterations", "Val");

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

    FilterType::Pointer iyFilter = FilterType::New();
    iyFilter->SetInput( petReader->GetOutput() );
    iyFilter->SetMaskInput( maskReader->GetOutput() );
    iyFilter->SetPSF(vVariance);
    iyFilter->SetIterations( nNumOfIters );
    iyFilter->SetVerbose ( bDebug );

    //Perform IY.
    try {
        iyFilter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tfailure applying discrete Iterative Yang on: " << sPETFileName
                  << "\n" << err
                  << std::endl;
        return EXIT_FAILURE;
    }

    PETWriterType::Pointer petWriter = PETWriterType::New();
    petWriter->SetFileName(sOutputFileName);
    petWriter->SetInput( iyFilter->GetOutput() );

    try {
        petWriter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "\n[Error]\tCannot write output file: " << sOutputFileName
                  << std::endl;

        return EXIT_FAILURE;
    }

    timer.Stop();

    if ( bDebug ) {
    	std::cout << "Completed in " << timer.GetTotal() << " seconds." << std::endl;
    }
    return EXIT_SUCCESS;
}

std::string getAcknowledgments(void)
{
    //Produces acknowledgments string for 3DSlicer.
    std::string sAck = "This program implements the Iterative Yang (IY) partial volume correction (PVC) technique. Please cite the following paper:\n"
                       "\tErlandsson, K. and Buvat, I. and Pretorius, P.H. and Thomas, B.A. and Hutton, B.F., (2012).\n\t\"A review of partial volume correction techniques "
                       "for emission tomography and their applications in neurology, cardiology and oncology\", \n\tPhysics in Medicine and Biology, vol. 57, no. 21, R119-59.";

    return sAck;
}

