#include <vector>
#include <string>

struct Bump
{
  double weight;
  double scale;
  double location[2];
};

typedef std::vector<Bump> LevelsetParameters;

double levelset(double x[], LevelsetParameters &params);

double heaviside(double x);

struct LevelsetDomain
{
  double xmin, xmax, ymin, ymax;
  double dx, dy;
};

void dump_levelset_parameters(LevelsetParameters &params);

void render_levelset_tofile(std::string filename, LevelsetParameters &params, LevelsetDomain grid);

