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
#include <cassert>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>

// command line parsing
#include "vul_arg.h"

// configuration file parsing
#include "json.h"

class RegisterParameterizedCFMApp
{
public:

  typedef short PixelType;
  typedef itk::Image<PixelType, 2> Image2DType;

  RegisterParameterizedCFMApp(int argc, char** argv)
    {
      // command line args
      vul_arg<std::string> target(0, "Target Input File");
      vul_arg<std::string> source(0, "Source Input File");
      vul_arg<std::string> output("-o", "Output File", "output.nii");
      vul_arg<std::string> conffile("-c", "Configuration File", "local.json");
      vul_arg_parse(argc, argv);

      targetFileName = target();
      sourceFileName = source();
      outputFileName = output();
      configurationFileName = conffile();

    }

  void parse_config()
    {
      Json::Value root;
      Json::Reader reader;
      bool parsingSuccessful = reader.parse( configurationFileName, root );
      if ( !parsingSuccessful )
	{
	std::cout  << "Failed to parse configuration\n"
		   << reader.getFormatedErrorMessages();
	}

      const Json::Value bumps = root["number-bumps"];
      numberOfBumpsX = bumps[Json::UInt(0)].asInt();
      numberOfBumpsY = bumps[Json::UInt(1)].asInt();

      initialSpatialScaleFactor = root["initial-spatial-scale-factor"].asDouble();
      initialIntensityScaleFactor = root["initial-intensity-scale-factor"].asDouble();
      regularizerWeight = root["regularizer-weight"].asDouble();
      matchingWeight = root["matching-weight"].asDouble();
      roisizeWeight = root["roisize-weight"].asDouble();
    }

  int run()
    {
      itk::ImageFileReader< Image2DType >::Pointer targetReader =
	itk::ImageFileReader< Image2DType >::New();
      targetReader->SetFileName(targetFileName);

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
    }

private:

  std::string targetFileName;
  std::string sourceFileName;
  std::string configurationFileName;
  std::string outputFileName;

  int numberOfBumpsX;
  int numberOfBumpsY;
  double initialSpatialScaleFactor;
  double initialIntensityScaleFactor;
  double regularizerWeight;
  double matchingWeight;
  double roisizeWeight;

};

int main(int argc, char** argv)
{
  RegisterParameterizedCFMApp app(argc, argv);

  app.parse_config();

  return app.run();

}
