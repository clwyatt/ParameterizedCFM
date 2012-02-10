#include <vector>
using std::vector;
#include <string>
using std::string;

#include "ParameterInitialization.h"
#include "ParameterizedLevelset.h"

const unsigned int DIMENSION = 2;

bool isodd(unsigned int x)
{
  return ((x % 2) == 1);
}

bool iseven(unsigned int x)
{
  return ((x % 2) == 0);
}

void initialize_levelset(ParameterOptions options, vnl_vector<double> &paramvector)
{
  unsigned int numberBumpsX = 2*options.numberCircles[0]+1;
  unsigned int numberBumpsY = 2*options.numberCircles[1]+1;

  double deltaX = (options.xmax - options.xmin)/((double)(numberBumpsX - 1));
  double deltaY = (options.ymax-options.ymin)/((double)(numberBumpsY-1));

  LevelsetParameters bumps;

  unsigned int bumpindex = 0;
  for(unsigned int row = 0; row < numberBumpsX; ++row)
    for(unsigned int col = 0; col < numberBumpsY; ++col)
      {
      if(iseven(row) && iseven(col)) continue;

      Bump b;
      bool positiveBump = isodd(row) && isodd(col);
      bool negativeBump = (iseven(row) && isodd(col))||(isodd(row) && iseven(col));

      double x = options.xmin + row*deltaX;
      double y = options.ymin + col*deltaY;

      if(positiveBump) b.weight = +options.weight;
      if(negativeBump) b.weight = -options.weight;
      b.scale = options.scale;
      b.location[0] = x;
      b.location[1] = y;

      bumps.push_back(b);
      bumpindex += 1;
      }


LevelsetDomain grid;
grid.xmin = options.xmin;
grid.xmax = options.xmax;
grid.ymin = options.ymin;
grid.ymax = options.ymax;
grid.dx = (options.xmax-options.xmin)/500.0;
grid.dy = (options.ymax-options.ymin)/500.0;

render_levelset_tofile(string("levelset.nrrd"), bumps, grid);

//rasterize into parametervector
unsigned int numberBumps = bumps.size();
unsigned int numberParameters = DIMENSION + numberBumps*(DIMENSION + 2);

  paramvector.set_size(numberParameters);
  paramvector[0] = 0; // x translation
  paramvector[1] = 0; // y translation
  // all alphas (weights)
  unsigned int paramidx = 2;
  for(unsigned int i = 0; i < numberBumps; ++i)
    {
    paramvector[paramidx] = bumps[i].weight;
    paramidx += 1;
    }
  // all betas (scales)
  for(unsigned int i = 0; i < numberBumps; ++i)
    {
    paramvector[paramidx] = bumps[i].scale;
    paramidx += 1;
    }
  // all chis (locations)
  for(unsigned int i = 0; i < numberBumps; ++i)
    {
    paramvector[paramidx] = bumps[i].location[0];
    paramvector[paramidx+1] = bumps[i].location[1];
    paramidx += 2;
    }
}
