/*****************************************************************************
Copyright (c) 2012, Bioimaging Systems Lab, Virginia Tech
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Virgina Tech nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*******************************************************************************/
#include <cstdlib>

#include <iostream>
using std::cerr;
using std::clog;
using std::endl;

// ITK
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include "itkImageFileWriter.h"

// CLI parsing
#include "vul_arg.h"

// define ITK short-hand types
typedef short PixelType;
typedef float InternalPixelType;
const unsigned int Dimension = 2;
typedef itk::Image< PixelType, Dimension >                              ImageType;
typedef itk::Image< InternalPixelType, Dimension >                      InternalImageType;
typedef itk::SmoothingRecursiveGaussianImageFilter< InternalImageType,
						    InternalImageType > FilterType;

const unsigned int IMAGE_SIZE = 64;
const double IMAGE_SPACING = 1.0;
const double CIRCLE_CENTER = 32;
const double CIRCLE_RADIUS = 24;
const double CIRCLE_INTENSITY = 4095.0;
const double HOLE_CENTER = 40;
const double HOLE_RADIUS = 5;
const double HOLE_INTENSITY = 0.0;
const double BACKGROUND_INTENSITY = 0.0;

// create a simple circle, possibly with another inside it
int makecircle(const char *filename, const float sigma,
	       const float shiftx, const float shifty, bool includehole)
{
  /// allocate image
  InternalImageType::Pointer   image = InternalImageType::New();
  InternalImageType::IndexType start;

  start[0] = 0; start[1] = 0;
  InternalImageType::SizeType size;
  size[0] = IMAGE_SIZE; size[1] = IMAGE_SIZE;
  InternalImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);
  image->SetRegions(region);
  InternalImageType::PointType origin;
  origin[0] = 0; origin[1] = 0;
  InternalImageType::SpacingType spacing;
  spacing[0] = IMAGE_SPACING; spacing[1] = IMAGE_SPACING;
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  typedef itk::ImageRegionIteratorWithIndex< InternalImageType > IteratorType;
  IteratorType regionIt( image, image->GetRequestedRegion() );
  for ( regionIt.GoToBegin(); !regionIt.IsAtEnd(); ++regionIt )
    {
    InternalImageType::IndexType idx = regionIt.GetIndex();
    if ( sqrt( pow(idx[0] - CIRCLE_CENTER - shiftx, 2) + pow(idx[1] - CIRCLE_CENTER - shifty, 2) ) < CIRCLE_RADIUS )
      {
      regionIt.Set(CIRCLE_INTENSITY);
      }
    else
      {
      regionIt.Set(BACKGROUND_INTENSITY);
      }
    if(includehole)
      {
      if ( sqrt( pow(idx[0] - HOLE_CENTER, 2) + pow(idx[1] - HOLE_CENTER, 2) ) < HOLE_RADIUS )
	{
	regionIt.Set(HOLE_INTENSITY);
	}
      }
    }

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->SetSigma(sigma);

  /// output the images
  typedef itk::ImageFileWriter< InternalImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput( filter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject & e )
    {
    cerr << "Output Image Write Failed." << endl;
    cerr << e << endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


int main(int argc, char *argv[])
{
  vul_arg< const char * > circleFilename(0, "Circle Image filename");
  vul_arg< const char * > circleWithHoleFilename(0, "Circle with Hole Image filename");
  vul_arg_parse(argc, argv);

  const float sigma = 1.0;

  if ( makecircle( circleFilename(), sigma, -5, -5, false ) == EXIT_FAILURE ) { return EXIT_FAILURE; }

  if ( makecircle( circleWithHoleFilename(), 0, 0, sigma, true ) == EXIT_FAILURE ) { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
