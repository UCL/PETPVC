/*
   IterativeYang.cxx
  
   Author:      Benjamin A. Thomas

   Copyright 2013 Institute of Nuclear Medicine, University College London.

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


#include <string>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkImageDuplicator.h>

#include <metaCommand.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_svd.h>
#include "FuzzyCorrFilter.h"

const char * const VERSION_NO="0.0.3";
const char * const AUTHOR="Benjamin A. Thomas";
const char * const APP_TITLE="Iterative Yang PVC";

typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<float, 4> MaskImageType;
typedef itk::Image<float, 3> PETImageType;

typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
typedef itk::ImageFileReader<PETImageType> PETReaderType;
typedef itk::ImageFileWriter<PETImageType> PETWriterType;


//Include fuzziness correction for the purposes of using soft masks.
typedef petpvc::FuzzyCorrFilter<MaskImageType> FuzzyCorrFilterType;

typedef itk::DiscreteGaussianImageFilter<PETImageType, PETImageType> BlurringFilterType;

//Extracts a 3D volume from 4D file.
typedef itk::ExtractImageFilter<MaskImageType, PETImageType> ExtractFilterType;

typedef itk::MultiplyImageFilter<PETImageType, PETImageType> MultiplyFilterType;
typedef itk::DivideImageFilter<PETImageType, PETImageType, PETImageType> DivideFilterType;
typedef itk::AddImageFilter<PETImageType, PETImageType> AddFilterType;

//For calculating mean values from image
typedef itk::StatisticsImageFilter<PETImageType> StatisticsFilterType;
typedef itk::ImageDuplicator<PETImageType> DuplicatorType;


//Function definitions:

//Creates image used to calculate correction factors.
PETImageType::Pointer getSyntheticPET(const MaskImageType::Pointer maskImage,
        const vnl_vector<float> vRegMeans);

//Generates the current estimate of the PET image for an iteration.
PETImageType::Pointer getCurrentEstimate(const PETImageType::Pointer origPET,
        const PETImageType::Pointer syntheticPET, BlurringFilterType::Pointer pBlurFilter);

//Produces the text for the acknowledgements dialog in Slicer. 
std::string getAcknowledgements(void);

int main(int argc, char *argv[]) {

    //Setting up command line argument list.
    MetaCommand command;

    command.SetVersion(VERSION_NO);
    command.SetAuthor(AUTHOR);
    command.SetName(APP_TITLE);
    command.SetDescription(
            "Performs iterative Yang (IY) partial volume correction");

    std::string sAcks = getAcknowledgements();
    command.SetAcknowledgments(sAcks.c_str());
    
    command.SetCategory("PETPVC");

    command.AddField("petfile", "PET filename", MetaCommand::IMAGE, MetaCommand::DATA_IN);
    command.AddField("maskfile", "mask filename", MetaCommand::FILE, MetaCommand::DATA_IN);
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
    command.AddOptionField("Iterations", "Val", MetaCommand::INT, false, "5");

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
    
    //Get number of iterations
    int nNumOfIters = command.GetValueAsInt("Iterations", "Val");

    //Make vector of FWHM in x,y and z.
    VectorType vFWHM;
    vFWHM[0] = fFWHM_x;
    vFWHM[1] = fFWHM_y;
    vFWHM[2] = fFWHM_z;

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

    //Create fuzziness correction filter.
    FuzzyCorrFilterType::Pointer fuzzyCorr = FuzzyCorrFilterType::New();
    fuzzyCorr->SetInput(maskReader->GetOutput());

    try {
        fuzzyCorr->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tFailed to calculate fuzziness correction factors!"
                << std::endl;
        return EXIT_FAILURE;
    }

    
    //Get fuzziness correction factors.
    vnl_matrix<float> fuzzCorrMat = fuzzyCorr->GetMatrix();
    vnl_vector<float> vecRegSize = fuzzyCorr->GetSumOfRegions();

    //Vector to contain the current estimate of the regional mean values.
    vnl_vector<float> vecRegMeansCurrent;
    vecRegMeansCurrent.set_size(vecRegSize.size());

    //Vector to contain the estimated means after fuzziness correction.
    vnl_vector<float> vecRegMeansUpdated;
    vecRegMeansUpdated.set_size(vecRegSize.size());

    //Calculate the variance for a given FWHM.
    VectorType vVariance;
    vVariance = vFWHM / (2.0 * sqrt(2.0 * log(2.0)));
    //std::cout << vVariance << std::endl;

    VectorType vVoxelSize = petReader->GetOutput()->GetSpacing();
    //std::cout << vVoxelSize << std::endl;

    //Scale for voxel size.
    vVariance[0] = pow((vVariance[0] / vVoxelSize[0]), 2);
    vVariance[1] = pow((vVariance[1] / vVoxelSize[1]), 2);
    vVariance[2] = pow((vVariance[2] / vVoxelSize[2]), 2);

    //Get mask image size.
    MaskImageType::SizeType imageSize =
            maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();

    int nClasses = 0;

    //If mask is not 4D, then quit.
    if (imageSize.Dimension == 4) {
        nClasses = imageSize[3];
    } else {
        std::cerr << "[Error]\tMask file: " << sMaskFileName << " must be 4-D!"
                << std::endl;
        return EXIT_FAILURE;
    }

    MaskImageType::IndexType desiredStart;
    desiredStart.Fill(0);
    MaskImageType::SizeType desiredSize = imageSize;

    //Create intermediate images to be used during iterations. 
    PETImageType::Pointer imageEstimate = PETImageType::New();
    PETImageType::Pointer imageSynthPET = PETImageType::New();
    PETImageType::Pointer imageExtractedRegion = PETImageType::New();

    //Extract filter used to extract 3D volume from 4D file.
    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput(maskReader->GetOutput());
    extractFilter->SetDirectionCollapseToIdentity(); // This is required.

    //Stats. filter used to calculate statistics for an image.
    StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
    
    //Multiplies two images together.
    MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    
    //Create blurring filter to apply PSF.
    BlurringFilterType::Pointer blurFilter = BlurringFilterType::New();
    blurFilter->SetVariance( vVariance );

    float fSumOfPETReg = 0.0;

    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(petReader->GetOutput());
    duplicator->Update();

    //Set image estimate to the original PET data for the first iteration.
    imageEstimate = duplicator->GetOutput();

    for (int k = 1; k <= nNumOfIters; k++) {
        if (k == 1) {
            std::cout << std::endl << "Iteration:  " << std::endl;
        }

        std::cout << k << "  "; // << std::endl;
        std::flush(std::cout);

        for (int i = 1; i <= nClasses; i++) {
            
            //Starts reading from 4D volume at index (0,0,0,i) through to 
            //(maxX, maxY, maxZ,0), i.e. one 3D brain mask.
            desiredStart[3] = i - 1;
            desiredSize[3] = 0;

            //Get region mask.
            extractFilter->SetExtractionRegion(
                    MaskImageType::RegionType(desiredStart, desiredSize));
            extractFilter->Update();

            imageExtractedRegion = extractFilter->GetOutput();
            imageExtractedRegion->SetDirection(imageEstimate->GetDirection());
            imageExtractedRegion->UpdateOutputData();

            
            //Multiply current image estimate by region mask. To clip PET values
            //to mask.
            multiplyFilter->SetInput1(imageEstimate);
            multiplyFilter->SetInput2(imageExtractedRegion);
            
            statsFilter->SetInput(multiplyFilter->GetOutput());
            statsFilter->Update();

            //Get sum of the clipped image.
            fSumOfPETReg = statsFilter->GetSum();
            
            //Place regional mean into vector.
            vecRegMeansCurrent.put(i - 1, fSumOfPETReg / vecRegSize.get(i - 1));
            //std::cout << "Sum = " << fSumOfPETReg << " , " << "Mean = " << vecRegMeansCurrent.get( i-1 ) << std::endl;

        }

        //Apply fuzziness correction to current mean value estimates.
        vecRegMeansUpdated = vnl_matrix_inverse<float>(fuzzCorrMat)
                * vecRegMeansCurrent;
        
        //vecRegMeansUpdated = vnl_svd< float >( fuzzCorrMat ).solve( vecRegMeansCurrent );
        //std::cout << vecRegMeansCurrent << std::endl;

        //Take the mask image and create pseudo PET image.
        imageSynthPET = getSyntheticPET(maskReader->GetOutput(),
                vecRegMeansUpdated);

        imageSynthPET->SetDirection(imageEstimate->GetDirection());
        imageSynthPET->UpdateOutputData();

        //Smooth pseudo image by PSF and take ratio: (pseudo image)/(smoothed).
        //The result of getCurrentEstimate() is the new image estimate for 
        //iteration k+1.
        imageEstimate = getCurrentEstimate(petReader->GetOutput(), imageSynthPET,
                        blurFilter);

    }

    std::cout << std::endl;

    //Write out result of final iteration.
    PETWriterType::Pointer petWriter = PETWriterType::New();
    petWriter->SetFileName(sOutputFileName);
    petWriter->SetInput(imageEstimate);

    try {
        petWriter->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "[Error]\tCannot write output file: " << sOutputFileName
                << std::endl;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

PETImageType::Pointer getSyntheticPET(const MaskImageType::Pointer maskImage,
        const vnl_vector<float> vRegMeans) {
    
    //Takes the 4D mask file along with the fuzziness-corrected mean values 
    //and creates the pseudo PET image.

    PETImageType::Pointer imageResult;// = PETImageType::New();

    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    AddFilterType::Pointer addFilter = AddFilterType::New();

    int nClasses = vRegMeans.size();

    MaskImageType::IndexType desiredStart;
    desiredStart.Fill(0);
    MaskImageType::SizeType desiredSize =
            maskImage->GetLargestPossibleRegion().GetSize();

    extractFilter->SetInput(maskImage);
    extractFilter->SetDirectionCollapseToIdentity();

    for (int i = 1; i <= nClasses; i++) {
        desiredStart[3] = i - 1;
        desiredSize[3] = 0;

        //Extract region mask.
        extractFilter->SetExtractionRegion(
                MaskImageType::RegionType(desiredStart, desiredSize));
        extractFilter->Update();

        //Multiply region mask by mean value.
        multiplyFilter->SetInput1(extractFilter->GetOutput());
        multiplyFilter->SetInput2(vRegMeans.get(i - 1));
        multiplyFilter->Update();

        //If this is the first region, create imageResult,
        //else add the current region to the previous contents of imageResult.
        if (i == 1) {
            imageResult = multiplyFilter->GetOutput();
            imageResult->DisconnectPipeline();
        } else {
            addFilter->SetInput1(imageResult);
            addFilter->SetInput2(multiplyFilter->GetOutput());
            addFilter->Update();

            imageResult = addFilter->GetOutput();
        }
    }

    return imageResult;
}

PETImageType::Pointer getCurrentEstimate(const PETImageType::Pointer origPET,
        const PETImageType::Pointer syntheticPET, BlurringFilterType::Pointer pBlurFilter) {

    //Takes the original PET data and the pseudo PET image, calculates the
    //correction factors  and returns the PV-corrected PET image.
    
    MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    DivideFilterType::Pointer divideFilter = DivideFilterType::New();

    //Smooth the pseudo PET by the PSF.
    //BlurringFilterType::Pointer blurringFilter = GaussianFilterType::New();
    pBlurFilter->SetInput(syntheticPET);
    //pBlurFilter->SetVariance(vVariance);

    //Take ratio of pseudo PET and smoothed pseudo PET. These are the correction
    //factors.
    divideFilter->SetInput1(syntheticPET);
    divideFilter->SetInput2( pBlurFilter->GetOutput() );
    
    //Multiply original PET by correction factors.
    multiplyFilter->SetInput1(origPET);
    multiplyFilter->SetInput2(divideFilter->GetOutput());
    multiplyFilter->Update();

    //Return new estimate.
    return multiplyFilter->GetOutput();

}

std::string getAcknowledgements(void) {
    //Produces acknowledgements string for 3DSlicer.
    std::string sAck = "This program implements the Iterative Yang (IY) partial volume correction (PVC) technique. Please cite the following paper:\n"
            "\tErlandsson, K. and Buvat, I. and Pretorius, P.H. and Thomas, B.A. and Hutton, B.F., (2012).\n\t\"A review of partial volume correction techniques "
            "for emission tomography and their applications in neurology, cardiology and oncology\", \n\tPhysics in Medicine and Biology, vol. 57, no. 21, R119-59.";

    return sAck;
}

