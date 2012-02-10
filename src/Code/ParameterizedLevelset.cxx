#include <iostream>
#include <cmath>

#include "ParameterizedLevelset.h"

#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkImageFileWriter.h"

const double EPS = 0.001;
const double PI = 3.14159265359;

double min(double a, double b)
{
  return (a < b) ? a : b;
}

double max(double a, double b)
{
  return (a < b) ? b : a;
}

template<unsigned int n>
double csrbf(double r)
{
  unsigned int el = (n >> 1) + 2;
  double temp = pow(max(0, 1-r), el+1);
  return temp*((el+1)*r+1);
}

double norm(double x[], double y[], double beta)
{
  double dx = x[0] - y[0];
  double dy = x[1] - y[1];
  return sqrt(beta*beta*dx*dx + beta*beta*dy*dy + EPS);
}

double levelset(double x[], LevelsetParameters &params)
{
  double sum = 0;
  for(unsigned int i = 0; i < params.size(); ++i)
    {
    Bump b = params[i];
    double alpha = b.weight;
    double beta = b.scale;

    double r = norm(x, b.location, beta);
    sum += (b.weight)*csrbf<2>(r);
    }
  return sum;
}


double heaviside(double x)
{
  return 0.5*(1 + (2/PI)*atan(PI*x/EPS));
}

void dump_levelset_parameters(LevelsetParameters &params)
{
  for(unsigned int i = 0; i < params.size(); ++i)
    {
    Bump b = params[i];
    std::cout << "[" << i << "] ";
    std::cout << b.weight << ", ";
    std::cout << b.scale << ", ";
    std::cout << b.location[0] << ", ";
    std::cout << b.location[1] << std::endl;
    }
}

void render_levelset_tofile(std::string filename, LevelsetParameters &params, LevelsetDomain grid)
{
  dump_levelset_parameters(params);

  typedef itk::Image<float, 2> InternalImageType;

  /// allocate image
  InternalImageType::Pointer   image = InternalImageType::New();
  InternalImageType::IndexType start;

  start[0] = 0; start[1] = 0;
  InternalImageType::SizeType size;
  size[0] = floor((grid.xmax-grid.xmin)/grid.dx);
  size[1] = floor((grid.ymax-grid.ymin)/grid.dy);
  InternalImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);
  image->SetRegions(region);
  InternalImageType::PointType origin;
  origin[0] = grid.xmin; origin[1] = grid.ymin;
  InternalImageType::SpacingType spacing;
  spacing[0] = grid.dx; spacing[1] = grid.dy;
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
  image->Allocate();

  typedef itk::ImageRegionIteratorWithIndex< InternalImageType > IteratorType;
  IteratorType regionIt( image, image->GetRequestedRegion() );
  double minimum = 1000000.0;
  double maximum = -1000000.0;
  for ( regionIt.GoToBegin(); !regionIt.IsAtEnd(); ++regionIt )
    {
    InternalImageType::IndexType idx = regionIt.GetIndex();
    double x[2];
    x[0] = origin[0] + spacing[0]*idx[0];
    x[1] = origin[1] + spacing[1]*idx[1];
    double value = heaviside(levelset(x, params));
    regionIt.Set(value);
    minimum = min(minimum, value);
    maximum = max(maximum, value);
    }

  std::cout << minimum << " to " << maximum << std::endl;

  /// output the image
  typedef itk::ImageFileWriter< InternalImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput( image );
  writer->Update();

}
