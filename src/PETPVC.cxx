/*
   PETPVC.cxx

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

 */

#include "EnvironmentInfo.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkCastImageFilter.h>
#include <metaCommand.h>
#include <vnl/vnl_matrix.h>

#include "petpvcRoussetPVCImageFilter.h"
#include "petpvcLabbePVCImageFilter.h"
#include "petpvcRBVPVCImageFilter.h"
#include "petpvcIterativeYangPVCImageFilter.h"
#include "petpvcMTCPVCImageFilter.h"
#include "petpvcMullerGartnerImageFilter.h"
#include "petpvcVanCittertPVCImageFilter.h"
#include "petpvcRLPVCImageFilter.h"
#include "petpvcSTCPVCImageFilter.h"

#include "petpvcLabbeRBVPVCImageFilter.h"
#include "petpvcLabbeMTCPVCImageFilter.h"

#include "petpvcIntraRegVCImageFilter.h"
#include "petpvcIntraRegRLImageFilter.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

enum PVCMethod { EGTM, ELabbe, EMullerGartner, EMTC,
				ERBV, EIterativeYang, ERichardsonLucy, EVanCittert,
				ELabbeRBV, ELabbeMTC, ERBVVanCittert, ERBVRichardsonLucy,
				EMTCVanCittert, EMTCRichardsonLucy, EIYVanCittert, EIYRichardsonLucy,
				ELabbeRBVVanCittert, ELabbeRBVRichardsonLucy, ELabbeMTCVanCittert, ELabbeMTCRichardsonLucy,
				EMGVanCittert, EMGRichardsonLucy, ESTC, EUnknown };

typedef itk::Vector<float, 3> VectorType;
typedef itk::Image<float, 4> MaskImageType;
typedef itk::Image<short, 3> Mask3DImageType;
typedef itk::Image<float, 3> PETImageType;

typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
typedef itk::ImageFileReader<Mask3DImageType> Mask3DReaderType;

typedef itk::ImageFileReader<PETImageType> PETReaderType;
typedef itk::ImageFileWriter<PETImageType> PETWriterType;

//Produces the text for the acknowledgment dialog in Slicer.
std::string getAcknowledgments(void);

//Takes input and checks for valid PVC method request.
PVCMethod getPVCMethod( std::string );

//Prints list of available methods.
void printPVCMethodList(void);

int main(int argc, char *argv[])
{

	const char * const AUTHOR = "Benjamin A. Thomas";
	const char * const APP_TITLE = "PETPVC";

    std::stringstream version_number;
    version_number << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
    const std::string VERSION_NO = version_number.str();

    //Setting up command line argument list.
    MetaCommand command;

    command.SetVersion(VERSION_NO.c_str());
    command.SetAuthor(AUTHOR);
    command.SetName(APP_TITLE);
    command.SetDescription(
        "PET partial volume correction toolbox");

    std::string sAcks = getAcknowledgments();
    command.SetAcknowledgments(sAcks.c_str());

    command.SetCategory("PETPVC");

    command.SetOption("Input", "i", true,"PET image file");
    command.SetOptionLongTag("Input", "input");
	command.AddOptionField("Input", "filename", MetaCommand::IMAGE, true, "");

    command.SetOption("Output", "o", true,"Output file");
    command.SetOptionLongTag("Output", "output");
	command.AddOptionField("Output", "filename", MetaCommand::STRING, true, "");

    command.SetOption("Mask", "m", false,"Mask image file");
    command.SetOptionLongTag("Mask", "mask");
	command.AddOptionField("Mask", "filename", MetaCommand::IMAGE, true, "");

    command.SetOption("PVC", "p", true,"Desired PVC method");
    command.SetOptionLongTag("PVC", "pvc");
	command.AddOptionField("PVC", "keyword", MetaCommand::STRING, true, "");

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

    command.SetOption("Iterations", "n", false, "Number of iterations");
    command.SetOptionLongTag("Iterations", "iter");
    command.AddOptionField("Iterations", "Val", MetaCommand::INT, false, "10");

    command.SetOption("Deconvolution", "k", false, "Number of deconvolution iterations");
    command.AddOptionField("Deconvolution", "Val", MetaCommand::INT, false, "10");

    command.SetOption("Alpha", "a", false, "Alpha value");
    command.SetOptionLongTag("Alpha", "alpha");
    command.AddOptionField("Alpha", "aval", MetaCommand::FLOAT, false, "1.5");

    command.SetOption("Stop", "s", false, "Stopping criterion");
    command.SetOptionLongTag("Stop", "stop");
    command.AddOptionField("Stop", "stopval", MetaCommand::FLOAT, false, "0.01");

    //Parse command line.
    if (!command.Parse(argc, argv)) {
		printPVCMethodList();
        return EXIT_FAILURE;
    }

    //Get image filenames
    std::string sPETFileName = command.GetValueAsString("Input", "filename");
    std::string sMaskFileName = command.GetValueAsString("Mask", "filename");
    std::string sOutputFileName = command.GetValueAsString("Output", "filename");

	std::string desiredMethod = command.GetValueAsString("PVC", "keyword");

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

	PVCMethod approach = getPVCMethod( desiredMethod );

	if (approach == EUnknown) {
		std::cerr << "[Error]\tUnknown method '" << desiredMethod << "' requested" << std::endl << std::endl;
		printPVCMethodList();
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

	//Create output image
	PETImageType::Pointer outputImage = NULL;

	//Create reader for mask image.
    MaskReaderType::Pointer maskReader = MaskReaderType::New();

	switch (approach) {
		case ERichardsonLucy: {
				std::cout << "Performing Richardson-Lucy..." << std::endl;

				typedef petpvc::RichardsonLucyPVCImageFilter< PETImageType >  RLFilterType;
				RLFilterType::Pointer rlFilter = RLFilterType::New();
    			rlFilter->SetInput( petReader->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfIters );

		    	//rlFilter->SetStoppingCond( fStop );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}
		case EVanCittert: {
				std::cout << "Performing reblurred Van-Cittert..." << std::endl;

				typedef petpvc::VanCittertPVCImageFilter< PETImageType >  VCFilterType;
				VCFilterType::Pointer vcFilter = VCFilterType::New();
    			vcFilter->SetInput( petReader->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();

				break;
			}
		default:
			maskReader->SetFileName(sMaskFileName);
			//Try to read mask.
    		try {
		        maskReader->Update();
		    } catch (itk::ExceptionObject & err) {
        		std::cerr << "[Error]\tCannot read mask input file: " << sMaskFileName
                  << std::endl;
        		return EXIT_FAILURE;
    		}
	}

    	switch (approach) {
		case ERBV: {
				std::cout << "Performing RBV..." << std::endl;
			    typedef petpvc::RBVPVCImageFilter<PETImageType, MaskImageType>  RBVFilterType;

				RBVFilterType::Pointer rbvFilter = RBVFilterType::New();
			    rbvFilter->SetInput( petReader->GetOutput() );
			    rbvFilter->SetMaskInput( maskReader->GetOutput() );
			    rbvFilter->SetPSF(vVariance);
			    rbvFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        rbvFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				outputImage = rbvFilter->GetOutput();

				break;
			}
		case EIterativeYang: {
				std::cout << "Performing iterative Yang..." << std::endl;
			    typedef petpvc::IterativeYangPVCImageFilter<PETImageType, MaskImageType>  IYFilterType;

				IYFilterType::Pointer iyFilter = IYFilterType::New();
			    iyFilter->SetInput( petReader->GetOutput() );
			    iyFilter->SetMaskInput( maskReader->GetOutput() );

				//Get number of iterations
				int nNumOfIters = command.GetValueAsInt("Iterations", "Val");
		    	iyFilter->SetIterations( nNumOfIters );

			    iyFilter->SetPSF(vVariance);
			    iyFilter->SetVerbose( bDebug );

			    //Perform iY.
			    try {
			        iyFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying iY on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				outputImage = iyFilter->GetOutput();

				break;
			}
		case EMTC: {
				std::cout << "Performing MTC..." << std::endl;
			    typedef petpvc::MTCPVCImageFilter<PETImageType, MaskImageType>  MTCFilterType;

				MTCFilterType::Pointer mtcFilter = MTCFilterType::New();
			    mtcFilter->SetInput( petReader->GetOutput() );
			    mtcFilter->SetMaskInput( maskReader->GetOutput() );
			    mtcFilter->SetPSF(vVariance);
			    mtcFilter->SetVerbose( bDebug );

			    //Perform MTC.
			    try {
			        mtcFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MTC on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				outputImage = mtcFilter->GetOutput();

				break;
			}
			case ESTC: {

					Mask3DReaderType::Pointer mask3Dreader = Mask3DReaderType::New();
					mask3Dreader->SetFileName(sMaskFileName);

					std::cout << "Performing STC..." << std::endl;
					typedef petpvc::STCPVCImageFilter<PETImageType, Mask3DImageType>  STCFilterType;

					STCFilterType::Pointer stcFilter = STCFilterType::New();
					stcFilter->SetInput( petReader->GetOutput() );
					stcFilter->SetMaskInput( mask3Dreader->GetOutput() );
					stcFilter->SetPSF(vVariance);
					stcFilter->SetVerbose( bDebug );

					//Perform STC.
					try {
						stcFilter->Update();
					} catch (itk::ExceptionObject & err) {
						std::cerr << "\n[Error]\tfailure applying STC on: " << sPETFileName
								  << "\n" << err
								  << std::endl;
						return EXIT_FAILURE;
					}

					outputImage = stcFilter->GetOutput();

					break;
		}
		case EMullerGartner: {
				std::cout << "Performing Muller-Gartner..." << std::endl;

				//Extracts a 3D volume from 4D file.
				typedef itk::ExtractImageFilter<MaskImageType, PETImageType> ExtractFilterType;

				MaskImageType::IndexType desiredStart;
    			desiredStart.Fill(0);
			    MaskImageType::SizeType desiredSize =
					maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();

			    //Extract filter used to extract 3D volume from 4D file.
			    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
			    extractFilter->SetInput(maskReader->GetOutput());
			    extractFilter->SetDirectionCollapseToIdentity(); // This is required.

			    //Extract filter used to extract 3D volume from 4D file.
			    ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
			    extractFilter2->SetInput(maskReader->GetOutput());
			    extractFilter2->SetDirectionCollapseToIdentity(); // This is required.

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

			    //Get WM mask.
				desiredStart[3] = 1;
			    extractFilter2->SetExtractionRegion(
			        MaskImageType::RegionType(desiredStart, desiredSize));
			    extractFilter2->Update();

			    imageWM = extractFilter2->GetOutput();
			    imageWM->SetDirection(petReader->GetOutput()->GetDirection());
			    imageWM->UpdateOutputData();

			    typedef petpvc::MullerGartnerImageFilter<PETImageType, PETImageType, PETImageType, PETImageType>  MGFilterType;

				MGFilterType::Pointer mgFilter = MGFilterType::New();
			    mgFilter->SetInput1(petReader->GetOutput());
    			mgFilter->SetInput2(imageGM);
    			mgFilter->SetInput3(imageWM);
				mgFilter->SetWM( 0 );
			    mgFilter->SetPSF(vVariance);
			    mgFilter->SetVerbose( bDebug );

			    //Perform MG.
			    try {
			        mgFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MG on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				outputImage = mgFilter->GetOutput();

				break;
			}

		case EMGVanCittert: {
				std::cout << "Performing Muller-Gartner..." << std::endl;

				//Extracts a 3D volume from 4D file.
				typedef itk::ExtractImageFilter<MaskImageType, PETImageType> ExtractFilterType;

				//Puts 3D into 4D.
				typedef itk::CastImageFilter<PETImageType, MaskImageType> CastFilterType;

				MaskImageType::IndexType desiredStart;
    			desiredStart.Fill(0);
			    MaskImageType::SizeType desiredSize =
					maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();

			    //Extract filter used to extract 3D volume from 4D file.
			    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
			    extractFilter->SetInput(maskReader->GetOutput());
			    extractFilter->SetDirectionCollapseToIdentity(); // This is required.

			    //Extract filter used to extract 3D volume from 4D file.
			    ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
			    extractFilter2->SetInput(maskReader->GetOutput());
			    extractFilter2->SetDirectionCollapseToIdentity(); // This is required.

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

			    //Get WM mask.
				desiredStart[3] = 1;
			    extractFilter2->SetExtractionRegion(
			        MaskImageType::RegionType(desiredStart, desiredSize));
			    extractFilter2->Update();

			    imageWM = extractFilter2->GetOutput();
			    imageWM->SetDirection(petReader->GetOutput()->GetDirection());
			    imageWM->UpdateOutputData();

			    typedef petpvc::MullerGartnerImageFilter<PETImageType, PETImageType, PETImageType, PETImageType>  MGFilterType;

				MGFilterType::Pointer mgFilter = MGFilterType::New();
			    mgFilter->SetInput1(petReader->GetOutput());
    			mgFilter->SetInput2(imageGM);
    			mgFilter->SetInput3(imageWM);
				mgFilter->SetWM( 0 );
			    mgFilter->SetPSF(vVariance);
			    mgFilter->SetVerbose( bDebug );

			    //Perform MG.
			    try {
			        mgFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MG on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Van-Cittert..." << std::endl;

				CastFilterType::Pointer castFilter = CastFilterType::New();
				castFilter->SetInput( imageGM );

				typedef petpvc::IntraRegVCImageFilter< PETImageType, MaskImageType >  IVCFilterType;
				IVCFilterType::Pointer vcFilter = IVCFilterType::New();
    			vcFilter->SetInput( mgFilter->GetOutput() );
				vcFilter->SetMaskInput( castFilter->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfDeconvIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();


				break;
			}
		case EMGRichardsonLucy: {
				std::cout << "Performing Muller-Gartner..." << std::endl;

				//Extracts a 3D volume from 4D file.
				typedef itk::ExtractImageFilter<MaskImageType, PETImageType> ExtractFilterType;

				//Puts 3D into 4D.
				typedef itk::CastImageFilter<PETImageType, MaskImageType> CastFilterType;

				MaskImageType::IndexType desiredStart;
    			desiredStart.Fill(0);
			    MaskImageType::SizeType desiredSize =
					maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();

			    //Extract filter used to extract 3D volume from 4D file.
			    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
			    extractFilter->SetInput(maskReader->GetOutput());
			    extractFilter->SetDirectionCollapseToIdentity(); // This is required.

			    //Extract filter used to extract 3D volume from 4D file.
			    ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
			    extractFilter2->SetInput(maskReader->GetOutput());
			    extractFilter2->SetDirectionCollapseToIdentity(); // This is required.

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

			    //Get WM mask.
				desiredStart[3] = 1;
			    extractFilter2->SetExtractionRegion(
			        MaskImageType::RegionType(desiredStart, desiredSize));
			    extractFilter2->Update();

			    imageWM = extractFilter2->GetOutput();
			    imageWM->SetDirection(petReader->GetOutput()->GetDirection());
			    imageWM->UpdateOutputData();

			    typedef petpvc::MullerGartnerImageFilter<PETImageType, PETImageType, PETImageType, PETImageType>  MGFilterType;

				MGFilterType::Pointer mgFilter = MGFilterType::New();
			    mgFilter->SetInput1(petReader->GetOutput());
    			mgFilter->SetInput2(imageGM);
    			mgFilter->SetInput3(imageWM);
				mgFilter->SetWM( 0 );
			    mgFilter->SetPSF(vVariance);
			    mgFilter->SetVerbose( bDebug );

			    //Perform MG.
			    try {
			        mgFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MG on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Richardson-Lucy..." << std::endl;

				CastFilterType::Pointer castFilter = CastFilterType::New();
				castFilter->SetInput( imageGM );

				typedef petpvc::IntraRegRLImageFilter< PETImageType, MaskImageType >  IRLFilterType;
				IRLFilterType::Pointer rlFilter = IRLFilterType::New();
    			rlFilter->SetInput( mgFilter->GetOutput() );
				rlFilter->SetMaskInput( castFilter->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfDeconvIters );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}
		case ELabbeRBV: {
				std::cout << "Performing Labbe-RBV..." << std::endl;
			    typedef petpvc::LabbeRBVPVCImageFilter<PETImageType, MaskImageType>  LabbeRBVFilterType;

				LabbeRBVFilterType::Pointer lrbvFilter = LabbeRBVFilterType::New();
			    lrbvFilter->SetInput( petReader->GetOutput() );
			    lrbvFilter->SetMaskInput( maskReader->GetOutput() );
			    lrbvFilter->SetPSF(vVariance);
			    lrbvFilter->SetVerbose( bDebug );

			    //Perform Labbe-RBV.
			    try {
			        lrbvFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying L-RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				outputImage = lrbvFilter->GetOutput();

				break;
			}
		case ELabbeMTC: {
				std::cout << "Performing Labbe-MTC..." << std::endl;
			    typedef petpvc::LabbeMTCPVCImageFilter<PETImageType, MaskImageType>  LabbeMTCFilterType;

				LabbeMTCFilterType::Pointer lmtcFilter = LabbeMTCFilterType::New();
			    lmtcFilter->SetInput( petReader->GetOutput() );
			    lmtcFilter->SetMaskInput( maskReader->GetOutput() );
			    lmtcFilter->SetPSF(vVariance);
			    lmtcFilter->SetVerbose( bDebug );

			    //Perform Labbe-MTC.
			    try {
			        lmtcFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying L-RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				outputImage = lmtcFilter->GetOutput();

				break;
			}

		case ERBVVanCittert: {
				std::cout << "Performing RBV..." << std::endl;
			    typedef petpvc::RBVPVCImageFilter<PETImageType, MaskImageType>  RBVFilterType;

				RBVFilterType::Pointer rbvFilter = RBVFilterType::New();
			    rbvFilter->SetInput( petReader->GetOutput() );
			    rbvFilter->SetMaskInput( maskReader->GetOutput() );
			    rbvFilter->SetPSF(vVariance);
			    rbvFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        rbvFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Van-Cittert..." << std::endl;

				typedef petpvc::IntraRegVCImageFilter< PETImageType, MaskImageType >  IVCFilterType;
				IVCFilterType::Pointer vcFilter = IVCFilterType::New();
    			vcFilter->SetInput( rbvFilter->GetOutput() );
				vcFilter->SetMaskInput( maskReader->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfDeconvIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();

				break;
			}
		case ERBVRichardsonLucy: {
				std::cout << "Performing RBV..." << std::endl;
			    typedef petpvc::RBVPVCImageFilter<PETImageType, MaskImageType>  RBVFilterType;

				RBVFilterType::Pointer rbvFilter = RBVFilterType::New();
			    rbvFilter->SetInput( petReader->GetOutput() );
			    rbvFilter->SetMaskInput( maskReader->GetOutput() );
			    rbvFilter->SetPSF(vVariance);
			    rbvFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        rbvFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Richardson-Lucy..." << std::endl;

				typedef petpvc::IntraRegRLImageFilter< PETImageType, MaskImageType >  IRLFilterType;
				IRLFilterType::Pointer rlFilter = IRLFilterType::New();
    			rlFilter->SetInput( rbvFilter->GetOutput() );
				rlFilter->SetMaskInput( maskReader->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfDeconvIters );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}

		case ELabbeRBVVanCittert: {
				std::cout << "Performing Labbe-RBV..." << std::endl;
			    typedef petpvc::LabbeRBVPVCImageFilter<PETImageType, MaskImageType>  LabbeRBVFilterType;

				LabbeRBVFilterType::Pointer lrbvFilter = LabbeRBVFilterType::New();
			    lrbvFilter->SetInput( petReader->GetOutput() );
			    lrbvFilter->SetMaskInput( maskReader->GetOutput() );
			    lrbvFilter->SetPSF(vVariance);
			    lrbvFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        lrbvFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Van-Cittert..." << std::endl;

				typedef petpvc::IntraRegVCImageFilter< PETImageType, MaskImageType >  IVCFilterType;
				IVCFilterType::Pointer vcFilter = IVCFilterType::New();
    			vcFilter->SetInput( lrbvFilter->GetOutput() );
				vcFilter->SetMaskInput( maskReader->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfDeconvIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();

				break;
			}
		case ELabbeRBVRichardsonLucy: {
				std::cout << "Performing Labbe-RBV..." << std::endl;
			    typedef petpvc::LabbeRBVPVCImageFilter<PETImageType, MaskImageType>  LabbeRBVFilterType;

				LabbeRBVFilterType::Pointer lrbvFilter = LabbeRBVFilterType::New();
			    lrbvFilter->SetInput( petReader->GetOutput() );
			    lrbvFilter->SetMaskInput( maskReader->GetOutput() );
			    lrbvFilter->SetPSF(vVariance);
			    lrbvFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        lrbvFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying RBV on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }


				std::cout << "Performing Intra-regional Richardson-Lucy..." << std::endl;

				typedef petpvc::IntraRegRLImageFilter< PETImageType, MaskImageType >  IRLFilterType;
				IRLFilterType::Pointer rlFilter = IRLFilterType::New();
    			rlFilter->SetInput( lrbvFilter->GetOutput() );
				rlFilter->SetMaskInput( maskReader->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfDeconvIters );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}

//MTC
		case EMTCVanCittert: {
				std::cout << "Performing MTC..." << std::endl;
			    typedef petpvc::MTCPVCImageFilter<PETImageType, MaskImageType>  MTCFilterType;

				MTCFilterType::Pointer mtcFilter = MTCFilterType::New();
			    mtcFilter->SetInput( petReader->GetOutput() );
			    mtcFilter->SetMaskInput( maskReader->GetOutput() );
			    mtcFilter->SetPSF(vVariance);
			    mtcFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        mtcFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MTC on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Van-Cittert..." << std::endl;

				typedef petpvc::IntraRegVCImageFilter< PETImageType, MaskImageType >  IVCFilterType;
				IVCFilterType::Pointer vcFilter = IVCFilterType::New();
    			vcFilter->SetInput( mtcFilter->GetOutput() );
				vcFilter->SetMaskInput( maskReader->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfDeconvIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();

				break;
			}
		case EMTCRichardsonLucy: {
				std::cout << "Performing MTC..." << std::endl;
			    typedef petpvc::MTCPVCImageFilter<PETImageType, MaskImageType>  MTCFilterType;

				MTCFilterType::Pointer mtcFilter = MTCFilterType::New();
			    mtcFilter->SetInput( petReader->GetOutput() );
			    mtcFilter->SetMaskInput( maskReader->GetOutput() );
			    mtcFilter->SetPSF(vVariance);
			    mtcFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        mtcFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MTC on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }


				std::cout << "Performing Intra-regional Richardson-Lucy..." << std::endl;

				typedef petpvc::IntraRegRLImageFilter< PETImageType, MaskImageType >  IRLFilterType;
				IRLFilterType::Pointer rlFilter = IRLFilterType::New();
    			rlFilter->SetInput( mtcFilter->GetOutput() );
				rlFilter->SetMaskInput( maskReader->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfDeconvIters );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}
		case ELabbeMTCVanCittert: {
				std::cout << "Performing Labbe-MTC..." << std::endl;
			    typedef petpvc::LabbeMTCPVCImageFilter<PETImageType, MaskImageType>  LabbeMTCFilterType;

				LabbeMTCFilterType::Pointer lmtcFilter = LabbeMTCFilterType::New();
			    lmtcFilter->SetInput( petReader->GetOutput() );
			    lmtcFilter->SetMaskInput( maskReader->GetOutput() );
			    lmtcFilter->SetPSF(vVariance);
			    lmtcFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        lmtcFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MTC on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Van-Cittert..." << std::endl;

				typedef petpvc::IntraRegVCImageFilter< PETImageType, MaskImageType >  IVCFilterType;
				IVCFilterType::Pointer vcFilter = IVCFilterType::New();
    			vcFilter->SetInput( lmtcFilter->GetOutput() );
				vcFilter->SetMaskInput( maskReader->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfDeconvIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();

				break;
			}
		case ELabbeMTCRichardsonLucy: {
				std::cout << "Performing Labbe-MTC..." << std::endl;
			    typedef petpvc::LabbeMTCPVCImageFilter<PETImageType, MaskImageType>  LabbeMTCFilterType;

				LabbeMTCFilterType::Pointer lmtcFilter = LabbeMTCFilterType::New();
			    lmtcFilter->SetInput( petReader->GetOutput() );
			    lmtcFilter->SetMaskInput( maskReader->GetOutput() );
			    lmtcFilter->SetPSF(vVariance);
			    lmtcFilter->SetVerbose( bDebug );

			    //Perform RBV.
			    try {
			        lmtcFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying MTC on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Richardson-Lucy..." << std::endl;

				typedef petpvc::IntraRegRLImageFilter< PETImageType, MaskImageType >  IRLFilterType;
				IRLFilterType::Pointer rlFilter = IRLFilterType::New();
    			rlFilter->SetInput( lmtcFilter->GetOutput() );
				rlFilter->SetMaskInput( maskReader->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfDeconvIters );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}
//IY
		case EIYVanCittert: {
				std::cout << "Performing iterative Yang..." << std::endl;
			    typedef petpvc::IterativeYangPVCImageFilter<PETImageType, MaskImageType>  IYFilterType;

				IYFilterType::Pointer iyFilter = IYFilterType::New();
			    iyFilter->SetInput( petReader->GetOutput() );
			    iyFilter->SetMaskInput( maskReader->GetOutput() );

				//Get number of iterations
				int nNumOfIters = command.GetValueAsInt("Iterations", "Val");
		    	iyFilter->SetIterations( nNumOfIters );

			    iyFilter->SetPSF(vVariance);
			    iyFilter->SetVerbose( bDebug );

			    //Perform iY.
			    try {
			        iyFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying iY on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Van-Cittert..." << std::endl;

				typedef petpvc::IntraRegVCImageFilter< PETImageType, MaskImageType >  IVCFilterType;
				IVCFilterType::Pointer vcFilter = IVCFilterType::New();
    			vcFilter->SetInput( iyFilter->GetOutput() );
				vcFilter->SetMaskInput( maskReader->GetOutput() );
		    	vcFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	vcFilter->SetIterations( nNumOfDeconvIters );

				//Get value for alpha.
    			float fAlpha = command.GetValueAsFloat("Alpha", "aval");
				vcFilter->SetAlpha( fAlpha );

    			//Get value for stopping criterion.
			    float fStop = command.GetValueAsFloat("Stop", "stopval");
		    	vcFilter->SetStoppingCond( fStop );

		    	vcFilter->SetVerbose ( bDebug );

    			//Perform VC.
    			try {
		    	    vcFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Van-Cittert on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = vcFilter->GetOutput();

				break;
			}
		case EIYRichardsonLucy: {
				std::cout << "Performing iterative Yang..." << std::endl;
			    typedef petpvc::IterativeYangPVCImageFilter<PETImageType, MaskImageType>  IYFilterType;

				IYFilterType::Pointer iyFilter = IYFilterType::New();
			    iyFilter->SetInput( petReader->GetOutput() );
			    iyFilter->SetMaskInput( maskReader->GetOutput() );

				//Get number of iterations
				int nNumOfIters = command.GetValueAsInt("Iterations", "Val");
		    	iyFilter->SetIterations( nNumOfIters );

			    iyFilter->SetPSF(vVariance);
			    iyFilter->SetVerbose( bDebug );

			    //Perform iY.
			    try {
			        iyFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying iY on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::cout << "Performing Intra-regional Richardson-Lucy..." << std::endl;

				typedef petpvc::IntraRegRLImageFilter< PETImageType, MaskImageType >  IRLFilterType;
				IRLFilterType::Pointer rlFilter = IRLFilterType::New();
    			rlFilter->SetInput( iyFilter->GetOutput() );
				rlFilter->SetMaskInput( maskReader->GetOutput() );
		    	rlFilter->SetPSF(vVariance);

				//Get number of iterations
				int nNumOfDeconvIters = command.GetValueAsInt("Deconvolution", "Val");
		    	rlFilter->SetIterations( nNumOfDeconvIters );
		    	rlFilter->SetVerbose ( bDebug );

    			//Perform RL.
    			try {
		    	    rlFilter->Update();
		    	} catch (itk::ExceptionObject & err) {
        			std::cerr << "[Error]\tfailure applying Richardson-Lucy on: " << sPETFileName
            	      << "\n" << err << std::endl;
        			return EXIT_FAILURE;
				}

				outputImage = rlFilter->GetOutput();

				break;
			}


			default:
				break;

		}

	if ( outputImage ) {
    	PETWriterType::Pointer petWriter = PETWriterType::New();
    	petWriter->SetFileName(sOutputFileName);
    	petWriter->SetInput( outputImage );

    	try {
    	    petWriter->Update();
    	} catch (itk::ExceptionObject & err) {
    	    std::cerr << "[Error]\tCannot write output file: " << sOutputFileName
                  << std::endl;

    	    return EXIT_FAILURE;
    	}
	}

	switch (approach) {
		case EGTM: {
				std::cout << "Performing Geometric matrix method..." << std::endl;
			    typedef petpvc::RoussetPVCImageFilter<PETImageType, MaskImageType>  GTMFilterType;

				GTMFilterType::Pointer gtmFilter = GTMFilterType::New();
			    gtmFilter->SetInput( petReader->GetOutput() );
			    gtmFilter->SetMaskInput( maskReader->GetOutput() );
			    gtmFilter->SetPSF(vVariance);
			    gtmFilter->SetVerbose( bDebug );

			    //Perform GTM.
			    try {
			        gtmFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying GTM on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::ofstream outputTextFile;
				outputTextFile.open( sOutputFileName.c_str() );
				if ( outputTextFile.is_open() ) {
					outputTextFile << "REGION\tMEAN" << std::endl;
					vnl_vector<float> results = gtmFilter->GetCorrectedMeans();
					for (int n=0; n < results.size(); n++)
						outputTextFile << n+1 << "\t" << results[n] << std::endl;
					outputTextFile.close();
				} else {
					 std::cerr << "[Error]\tCannot write output file: " << sOutputFileName
                  		<< std::endl;

    	    		return EXIT_FAILURE;
				}

				break;
			}
		case ELabbe: {
				std::cout << "Performing the Labbe method..." << std::endl;
			    typedef petpvc::LabbePVCImageFilter<PETImageType, MaskImageType>  LabbeFilterType;

				LabbeFilterType::Pointer labbeFilter = LabbeFilterType::New();
			    labbeFilter->SetInput( petReader->GetOutput() );
			    labbeFilter->SetMaskInput( maskReader->GetOutput() );
			    labbeFilter->SetPSF(vVariance);
			    labbeFilter->SetVerbose( bDebug );

			    //Perform L-PVC.
			    try {
			        labbeFilter->Update();
			    } catch (itk::ExceptionObject & err) {
			        std::cerr << "\n[Error]\tfailure applying Labbe on: " << sPETFileName
			                  << "\n" << err
			                  << std::endl;
			        return EXIT_FAILURE;
			    }

				std::ofstream outputTextFile;
				outputTextFile.open( sOutputFileName.c_str() );
				if ( outputTextFile.is_open() ) {
					outputTextFile << "REGION\tMEAN" << std::endl;
					vnl_vector<float> results = labbeFilter->GetCorrectedMeans();
					for (int n=0; n < results.size(); n++)
						outputTextFile << n+1 << "\t" << results[n] << std::endl;
					outputTextFile.close();
				} else {
					 std::cerr << "[Error]\tCannot write output file: " << sOutputFileName
                  		<< std::endl;

    	    		return EXIT_FAILURE;
				}

				break;
			}
		default: break;
	}


    return EXIT_SUCCESS;
}

std::string getAcknowledgments(void)
{
    //Produces acknowledgments string for 3DSlicer.
    std::string sAck = "";
    return sAck;
}


PVCMethod getPVCMethod( std::string method ) {

	std::transform(method.begin(), method.end(),method.begin(), ::toupper);

	if ( method == "GTM" )
		return EGTM;

	if ( method == "LABBE" )
		return ELabbe;

	if ( method == "MG" )
		return EMullerGartner;

	if ( method == "MTC" )
		return EMTC;

	if ( method == "RBV" )
		return ERBV;

	if ( method == "IY" )
		return EIterativeYang;

	if ( method == "RL" )
		return ERichardsonLucy;

	if ( method == "STC" )
		return ESTC;

	if ( method == "VC" )
		return EVanCittert;

	if ( ( method == "LABBE+RBV" ) || ( method == "RBV+LABBE" ) )
		return ELabbeRBV;

	if ( ( method == "LABBE+MTC" ) || ( method == "MTC+LABBE" ) )
		return ELabbeMTC;

	if ( method == "RBV+VC" )
		return ERBVVanCittert;

	if ( method == "RBV+RL" )
		return ERBVRichardsonLucy;

	if ( method == "LABBE+RBV+VC" )
		return ELabbeRBVVanCittert;

	if ( method == "LABBE+RBV+RL" )
		return ELabbeRBVRichardsonLucy;

	if ( method == "MTC+VC" )
		return EMTCVanCittert;

	if ( method == "MTC+RL" )
		return EMTCRichardsonLucy;

	if ( method == "LABBE+MTC+VC" )
		return ELabbeMTCVanCittert;

	if ( method == "LABBE+MTC+RL" )
		return ELabbeMTCRichardsonLucy;

	if ( method == "IY+VC" )
		return EIYVanCittert;

	if ( method == "IY+RL" )
		return EIYRichardsonLucy;

	if ( method == "MG+VC" )
		return EMGVanCittert;

	if ( method == "MG+RL" )
		return EMGRichardsonLucy;

	return EUnknown;
}

void printPVCMethodList(void) {

	std::cout << std::endl << "----------------------------------------------" << std::endl;
	std::cout << "Technique - keyword" << std::endl << std::endl;

	std::cout << "Geometric transfer matrix - \"GTM\"" << std::endl;
	std::cout << "Labbe approach - \"LABBE\"" << std::endl;
	std::cout << "Richardson-Lucy - \"RL\"" << std::endl;
	std::cout << "Van-Cittert - \"VC\"" << std::endl;

	std::cout << "Region-based voxel-wise correction - \"RBV\"" << std::endl;
	std::cout << "RBV with Labbe - \"LABBE+RBV\"" << std::endl;
	std::cout << "RBV with Van-Cittert - \"RBV+VC\"" << std::endl;
	std::cout << "RBV with Richardson-Lucy - \"RBV+RL\"" << std::endl;
	std::cout << "RBV with Labbe and Van-Cittert - \"LABBE+RBV+VC\"" << std::endl;
	std::cout << "RBV with Labbe and Richardson-Lucy- \"LABBE+RBV+RL\"" << std::endl;

	std::cout << "Single-target correction - \"STC\"" << std::endl;

	std::cout << "Multi-target correction - \"MTC\"" << std::endl;
	std::cout << "MTC with Labbe - \"LABBE+MTC\"" << std::endl;
	std::cout << "MTC with Van-Cittert - \"MTC+VC\"" << std::endl;
	std::cout << "MTC with Richardson-Lucy - \"MTC+RL\"" << std::endl;
	std::cout << "MTC with Labbe and Van-Cittert - \"LABBE+MTC+VC\"" << std::endl;
	std::cout << "MTC with Labbe and Richardson-Lucy- \"LABBE+MTC+RL\"" << std::endl;

	std::cout << "Iterative Yang - \"IY\"" << std::endl;
	std::cout << "Iterative Yang with Van-Cittert - \"IY+VC\"" << std::endl;
	std::cout << "Iterative Yang with Richardson-Lucy - \"IY+RL\"" << std::endl;

	std::cout << "Muller Gartner - \"MG\"" << std::endl;
	std::cout << "Muller Gartner with Van-Cittert - \"MG+VC\"" << std::endl;
	std::cout << "Muller Gartner with Richardson-Lucy - \"MG+RL\"" << std::endl;


	std::cout << std::endl;
}
