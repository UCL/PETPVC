/*=========================================================================
 *  Author:      Kris Thielemans
 *
 *  Copyright Insight Software Consortium
 *  Copyright University College London
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/*
   This file creates a test image consisting of several (uniform) geometric
   objects, and its corresponding segmentation.

   The code to create the image is largely based on ITK's
      Examples/Filtering/SpatialObjectToImage1.cxx

   Update: Nov. 08 2014 (V. Cuplov)
   In order to test Ben's new 'discrete' Iterative Yang code pvc_diy, which only 
   needs a 3D image with labels (a parcellation image) as opposed to the 4D mask, 
   this file now creates an additionnal image: 3D parcellation of the original image 
   in which the regions are labelled with some random values.

*/

#include "itkSpatialObjectToImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkCylinderSpatialObject.h"
#include "itkGroupSpatialObject.h"
#include "itkImageFileWriter.h"
#include <itkBinaryThresholdImageFilter.h>

#include "itkCastImageFilter.h"
#include "itkPasteImageFilter.h"
#include <iostream>

// helper function:
// store a 3D image inside a 4D image at a particular index
// at present, we return a new 4D image (as KT couldn't figure out how to run "in-place")
template <typename ImageType, typename Image4DType>
typename Image4DType::Pointer
copy_3d_to_4d(const ImageType * image, const Image4DType * image4d, int idx)
{
    // this works by first casting the 3d image to a 4d one (with a single "frame")
    // then "pasting" this 4d image inside the image4d at the correct place

    typedef itk::CastImageFilter<ImageType, Image4DType> castType;
    typedef itk::PasteImageFilter<Image4DType, Image4DType> pasteType;

    typename castType::Pointer castFilter = castType::New();
    castFilter->SetInput(image);
    castFilter->Update();
    typename pasteType::Pointer pasteFilter = pasteType::New();
    pasteFilter->SetSourceImage(castFilter->GetOutput());
    pasteFilter->SetDestinationImage(image4d);
    pasteFilter->SetSourceRegion(castFilter->GetOutput()->GetLargestPossibleRegion());
    //pasteFilter->InPlaceOn(); // this doesn't work for some reason
    typename Image4DType::IndexType destinationIndex;
    destinationIndex[0] = 0;
    destinationIndex[1] = 0;
    destinationIndex[2] = 0;
    destinationIndex[3] = idx;
    pasteFilter->SetDestinationIndex(destinationIndex);
    pasteFilter->Update();
    typename Image4DType::Pointer newimage4d = pasteFilter->GetOutput();
    return newimage4d;
}

int main( int argc, char *argv[] )
{
    if( argc != 4 ) {
        std::cerr << "Usage: " << argv[0] << " outputimagefile outputmaskfile outputparcellation" << std::endl;
        return EXIT_FAILURE;
    }
    //  We declare the pixel type and dimension of the image to be produced as
    //  output.
    typedef signed short  PixelType;
    const unsigned int    Dimension = 3;

    typedef itk::Image< PixelType, Dimension >       ImageType;
    typedef itk::Image< PixelType, Dimension+1 >     MaskType;
    // size of the output image
    // (let's make it asymmetric for some additional testing...)
    ImageType::SizeType size;
    size[ 0 ] =  50;
    size[ 1 ] =  40;
    size[ 2 ] = 60;
    ImageType::SpacingType spacing;
    spacing[0] =  100.0 / size[0];
    spacing[1] =  100.0 / size[1];
    spacing[2] =  300.0 / size[2];

    try {
        //  we instantiate the types of the elementary
        //  SpatialObjects that we plan to group, and we instantiate as well the type
        //  of the SpatialObject that will hold the group together.
        typedef itk::EllipseSpatialObject< Dimension >   EllipseType;
        typedef itk::CylinderSpatialObject               CylinderType;
        typedef itk::GroupSpatialObject< Dimension >     GroupType;
        typedef itk::SpatialObject< Dimension >          SpatialObjectType;
        //
        //  We instantiate the SpatialObjectToImageFilter type by using as template
        //  arguments the input SpatialObject and the output image types.
        //
        typedef itk::SpatialObjectToImageFilter<
        SpatialObjectType, ImageType >   SpatialObjectToImageFilterType;

        SpatialObjectToImageFilterType::Pointer imageFilter =
            SpatialObjectToImageFilterType::New();

        imageFilter->SetSize( size );
        imageFilter->SetSpacing( spacing );

        // we will store the original image in this variable
        ImageType::Pointer image;
        // we will store the parcellation image in this variable
        ImageType::Pointer parcellation;

        //  We create the elementary shapes that are going to be composed into the
        //  group spatial objects.
        EllipseType::Pointer ellipse    = EllipseType::New();
        CylinderType::Pointer cylinder1 = CylinderType::New();
        CylinderType::Pointer cylinder2 = CylinderType::New();
        const PixelType backgroundValue  = 0;
        const PixelType ellipseValue = 2;
        const PixelType cylinder1Value = 3;
        const PixelType cylinder2Value = 4;

// for parcellation image
        const PixelType ellipseparcellationValue = 124;
        const PixelType cylinder1parcellationValue = 17;
        const PixelType cylinder2parcellationValue = 49;


        // create these objects and the corresponding image
        {
            //
            //  The Elementary shapes have internal parameters of their own. These
            //  parameters define the geometrical characteristics of the basic shapes.
            //  For example, a cylinder is defined by its radius and height.
            ellipse->SetRadius(  size[0] * 0.2 * spacing[0] );

            cylinder1->SetRadius(  size[0] * 0.2 * spacing[0] );
            cylinder2->SetRadius(  size[0] * 0.2 * spacing[0] );

            cylinder1->SetHeight( size[2] * 0.30 * spacing[2]);
            cylinder2->SetHeight( size[2] * 0.30 * spacing[2]);
            //
            //  Each one of these components will be placed in a different position and
            //  orientation. We define transforms in order to specify those relative
            //  positions and orientations.
            typedef GroupType::TransformType                 TransformType;

            TransformType::Pointer transform1 = TransformType::New();
            TransformType::Pointer transform2 = TransformType::New();
            TransformType::Pointer transform3 = TransformType::New();

            transform1->SetIdentity();
            transform2->SetIdentity();
            transform3->SetIdentity();
            //
            //  Then we set the specific values of the transform parameters, and we
            //  assign the transforms to the elementary shapes.
            TransformType::OutputVectorType  translation;
            TransformType::CenterType        center;

            translation[ 0 ] =  size[0] * spacing[0] / 2.0;
            translation[ 1 ] =  size[1] * spacing[1] / 4.0;
            translation[ 2 ] =  size[2] * spacing[2] / 2.0;
            transform1->Translate( translation, false );

            translation[ 1 ] =  size[1] * spacing[1] / 2.0;
            translation[ 2 ] =  size[2] * spacing[2] * 0.22;
            transform2->Rotate( 1, 2, vnl_math::pi / 2.0 );
            transform2->Translate( translation, false );

            translation[ 2 ] = size[2] * spacing[2] * 0.78;
            transform3->Rotate( 1, 2, vnl_math::pi / 2.0 );
            transform3->Translate( translation, false );

            ellipse->SetObjectToParentTransform( transform1 );
            cylinder1->SetObjectToParentTransform( transform2 );
            cylinder2->SetObjectToParentTransform( transform3 );
            //
            //  The elementary shapes are aggregated in a parent group, that in turn is
            //  passed as input to the filter.
            GroupType::Pointer group = GroupType::New();
            group->AddSpatialObject( ellipse );
            group->AddSpatialObject( cylinder1 );
            group->AddSpatialObject( cylinder2 );

            //
            //  By default, the filter will rasterize the aggregation of elementary
            //  shapes and will assign a pixel value to locations that fall inside of any
            //  of the elementary shapes, and a different pixel value to locations that
            //  fall outside of all of the elementary shapes. It is possible, however, to
            //  generate richer images if we allow the filter to use the values that the
            //  elementary spatial objects return via their \code{ValueAt} methods. This
            //  is what we choose to do in this example, by using the following code.

            ellipse->SetDefaultInsideValue( ellipseValue   );
            cylinder1->SetDefaultInsideValue( cylinder1Value );
            cylinder2->SetDefaultInsideValue( cylinder2Value );

            ellipse->SetDefaultOutsideValue(   backgroundValue );
            cylinder1->SetDefaultOutsideValue( backgroundValue );
            cylinder2->SetDefaultOutsideValue( backgroundValue );

            imageFilter->SetInput(  group  );
            imageFilter->SetUseObjectValue( true );
            imageFilter->SetOutsideValue( backgroundValue );
            imageFilter->Update();
            // keep image for later
            image = imageFilter->GetOutput();

        // write to file
            typedef itk::ImageFileWriter< ImageType >     WriterType;
            WriterType::Pointer writer = WriterType::New();

            writer->SetFileName( argv[1] );
            writer->SetInput( image);

            writer->Update();

// values for regions in the parcellation image:
            ellipse->SetDefaultInsideValue( ellipseparcellationValue   );
            cylinder1->SetDefaultInsideValue( cylinder1parcellationValue );
            cylinder2->SetDefaultInsideValue( cylinder2parcellationValue );

            imageFilter->SetInput(  group  );
            imageFilter->SetUseObjectValue( true );
            imageFilter->SetOutsideValue( backgroundValue );
            imageFilter->Update();

            // keep parcellation image for later
            parcellation = imageFilter->GetOutput();

        // write to file
            writer->SetFileName( argv[3] );
            writer->SetInput( parcellation );
            writer->Update();

        }


        // now create mask: it has to be a 4D image (each volume corresponding to a region)
        MaskType::Pointer mask  = MaskType::New();
        // create a 4D image with appropriate geometry
        {
            MaskType::IndexType start;
            for (int i=0; i<3; ++i) start[i] = 0;
            start[3] = 0;
            MaskType::SizeType masksize;
            for (int i=0; i<3; ++i) masksize[i] = size[i];
            masksize[3] = 4; // background + 3 regions
            MaskType::RegionType region;
            region.SetSize( masksize );
            region.SetIndex( start );
            mask->SetRegions( region );

            // tediously copy spatial info
            // (there might be a better way...)
            MaskType::DirectionType direction;
            for (int i=0; i<3; ++i) {
                for (int j=0; j<3; ++j)
                    direction[i][j]=image->GetDirection()[i][j];
                direction[3][i]=0;
                direction[i][3]=0;
            }
            direction[3][3]=1;
            mask->SetDirection(direction);
            MaskType::SpacingType spacing;
            for (int i=0; i<3; ++i) spacing[i]=image->GetSpacing()[i];
            spacing[3]=1;
            mask->SetSpacing(spacing);
            MaskType::PointType origin;
            for (int i=0; i<3; ++i) origin[i]=image->GetOrigin()[i];
            origin[3]=0;
            mask->SetOrigin(origin);
            mask->Allocate();
        }

        // now fill the mask appropriately
        {
            // first find background, store that as first mask
            {
                typedef itk::BinaryThresholdImageFilter <ImageType, ImageType> BinThresholdImageFilterType;
                BinThresholdImageFilterType::Pointer thresholdFilter = BinThresholdImageFilterType::New();
                thresholdFilter->SetInput( image );
                thresholdFilter->SetUpperThreshold( backgroundValue );
                thresholdFilter->SetLowerThreshold( backgroundValue );
                thresholdFilter->SetInsideValue( 1 );
                thresholdFilter->SetOutsideValue( 0 );
                thresholdFilter->Update();
                // copy to mask
                mask=copy_3d_to_4d(thresholdFilter->GetOutput(), mask.GetPointer(), 0);
            }

            // find the mask for each object, just by reusing the original object as
            // input for imageFilter
            // (we could have used thresholding as well, as long as all values are different)
            imageFilter->SetUseObjectValue( false );
            imageFilter->SetOutsideValue( 0 );
            imageFilter->SetInsideValue( 1 ); // "in the object" has to be 1

            imageFilter->SetInput( ellipse );
            imageFilter->Update();
            mask=copy_3d_to_4d(imageFilter->GetOutput(), mask.GetPointer(), 1);

            imageFilter->SetInput( cylinder1 );
            imageFilter->Update();
            mask=copy_3d_to_4d(imageFilter->GetOutput(), mask.GetPointer(), 2);

            imageFilter->SetInput( cylinder2 );
            imageFilter->Update();
            mask=copy_3d_to_4d(imageFilter->GetOutput(), mask.GetPointer(), 3);
        }


        // now write it do disk!
        {
            typedef itk::ImageFileWriter< MaskType >     MaskWriterType;
            MaskWriterType::Pointer maskWriter = MaskWriterType::New();
            maskWriter->SetFileName( argv[2] );
            maskWriter->SetInput( mask );
            maskWriter->Update();
        }
    } catch( itk::ExceptionObject & excp ) {
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
