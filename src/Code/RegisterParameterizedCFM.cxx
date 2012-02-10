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
#include <iostream>
using std::cout;
using std::cerr;
using std::clog;
using std::endl;

#include <string>
using std::string;

#include <fstream>
using std::ifstream;

// command line parsing
#include "vul_arg.h"

// configuration file parsing
#include "json.h"

// numerics
#include <vnl/vnl_vector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>

#include "ParameterInitialization.h"

class RegisterParameterizedCFMApp
{
public:

  typedef short PixelType;
  typedef itk::Image<PixelType, 2> Image2DType;

  RegisterParameterizedCFMApp(int argc, char** argv)
    {
      // command line args
      vul_arg<string> target(0, "Target Input File");
      vul_arg<string> source(0, "Source Input File");
      vul_arg<string> output("-o", "Output File", "output.nii");
      vul_arg<string> conffile("-c", "Configuration File", "local.json");
      vul_arg_parse(argc, argv);

      targetFileName = target();
      sourceFileName = source();
      outputFileName = output();
      configurationFileName = conffile();

    }

  // TODO this function is tooooo loooooong
  bool parse_config()
    {
      ifstream infile(configurationFileName.c_str());
      if( infile.fail() )
	{
	cout << "Cannot open configuration file." << endl;
	return false;
	}

      Json::Value root;
      Json::Reader reader;
      bool parsingSuccessful = reader.parse( infile, root );
      infile.close();

      if ( !parsingSuccessful )
	{
	std::cout  << "Failed to parse configuration\n"
		   << reader.getFormatedErrorMessages();
	return false;
	}

      // TODO add error handling below
      if(root.isMember("number-bumps"))
	{
	Json::UInt index(0);
	if(!root["number-bumps"].isValidIndex(index) )
	  {
	  cout << "Warning: number-bumps-x not in configuration." << endl;
	  }
	index = 1;
	if(!root["number-bumps"].isValidIndex(index) )
	  {
	  cout << "Warning: number-bumps-y not in configuration." << endl;
	  }
	}
      else
	{
	cout << "Warning: number-bumps not in configuration." << endl;
	}

      if(!root.isMember("initial-spatial-scale-factor"))
	{
	cout << "Warning: initial-spatial-scale-factor not in configuration." << endl;
	}
      if(!root.isMember("initial-intensity-scale-factor"))
	{
	cout << "Warning: initial-intensity-scale-factor not in configuration." << endl;
	}
      if(!root.isMember("regularizer-weight"))
	{
	cout << "Warning: regularizer-weight not in configuration." << endl;
	}
      if(!root.isMember("matching-weight"))
	{
	cout << "Warning: matching-weight not in configuration." << endl;
	}
      if(!root.isMember("roisize-weight"))
	{
	cout << "Warning: roisize-weight not in configuration." << endl;
	}

      // read all using defaults if not present
      const Json::Value bumps = root["number-bumps"];
      numberOfBumpsX = bumps[Json::UInt(0)].asInt();
      numberOfBumpsY = bumps[Json::UInt(1)].asInt();
      initialSpatialScaleFactor = root["initial-spatial-scale-factor"].asDouble();
      initialIntensityScaleFactor = root["initial-intensity-scale-factor"].asDouble();
      regularizerWeight = root["regularizer-weight"].asDouble();
      matchingWeight = root["matching-weight"].asDouble();
      roisizeWeight = root["roisize-weight"].asDouble();

      return true;
    }

  int run()
    {
      itk::ImageFileReader< Image2DType >::Pointer targetReader =
	itk::ImageFileReader< Image2DType >::New();
      targetReader->SetFileName(targetFileName);

      try
	{
	targetReader->Update();
	target = targetReader->GetOutput();
	}
      catch( itk::ExceptionObject & ex )
	{
	std::cerr << "Exception thrown while reading target.\n";
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

      computeInitialParameters();
      cout << "Number of optimization parameters = " << initialParameters.size() << endl;

      itk::ImageFileReader< Image2DType >::Pointer sourceReader =
	itk::ImageFileReader< Image2DType >::New();
      sourceReader->SetFileName(sourceFileName);

      itk::ImageFileWriter< Image2DType >::Pointer writer =
	itk::ImageFileWriter< Image2DType >::New();
      writer->SetFileName(outputFileName);

      writer->SetInput( sourceReader->GetOutput() );

      try
	{
	writer->Update();
	}
      catch( itk::ExceptionObject & excp )
	{
	std::cerr << "Exception thrown while writing.\n";
	std::cerr << excp << std::endl;
	return EXIT_FAILURE;
	}

      return EXIT_SUCCESS;
    }

  void computeInitialParameters()
    {
      // get target image extents
      Image2DType::RegionType targetRegion = target->GetLargestPossibleRegion();
      Image2DType::SizeType targetSize = targetRegion.GetSize();
      Image2DType::SpacingType targetSpacing = target->GetSpacing();
      Image2DType::PointType targetOrigin = target->GetOrigin();
      double xmin = targetOrigin[0];
      double xmax = xmin+targetSize[0]*targetSpacing[0];
      double ymin = targetOrigin[1];
      double ymax = ymin+targetSize[1]*targetSpacing[1];

      // compute initial parameters
      ParameterOptions options;
      options.numberCircles[0] = numberOfBumpsX;
      options.numberCircles[1] = numberOfBumpsY;
      options.weight = initialIntensityScaleFactor;
      options.scale = initialSpatialScaleFactor;
      options.xmin = xmin;
      options.xmax = xmax;
      options.ymin = ymin;
      options.ymax = ymax;

      initialize_levelset(options, initialParameters);

    }

private:

  std::string targetFileName;
  std::string sourceFileName;
  std::string configurationFileName;
  std::string outputFileName;

  Image2DType::Pointer target;

  int numberOfBumpsX;
  int numberOfBumpsY;
  double initialSpatialScaleFactor;
  double initialIntensityScaleFactor;
  double regularizerWeight;
  double matchingWeight;
  double roisizeWeight;

  vnl_vector<double> initialParameters;

};

int main(int argc, char** argv)
{
  RegisterParameterizedCFMApp app(argc, argv);

  if(!app.parse_config()) return EXIT_FAILURE;

  return app.run();

}
