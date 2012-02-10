#include <vnl/vnl_vector.h>

struct ParameterOptions
{
  double numberCircles[2];
  double weight;
  double scale;
  double xmin, xmax, ymin, ymax;
};

void initialize_levelset(ParameterOptions options, vnl_vector<double> &paramvector);
