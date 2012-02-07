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

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <string>
using std::string;

#include <cassert>

#include <vcl_iosfwd.h>
#include <vnl/algo/vnl_lbfgs.h>
#include <vnl/vnl_vector.h>

// for command line parsing
#include "vul_arg.h"

typedef double SampleType;
typedef vnl_vector<SampleType> SignalType;

const double EPS = 0.001;
const double PI = 3.14159265359;

void linspace(SignalType & x, double &step, double low, double high, unsigned int sz)
{
  double val = low;
  step = (high - low)/(static_cast<double>(sz)-1);
  x.set_size(sz);
  for(unsigned int i = 0; i < x.size(); ++i)
    {
    x[i] = val;
    val += step;
    }
}

unsigned int locate(SignalType & xx, double x)
{
  unsigned int low, up;
  low = 0;
  up = xx.size() - 1;
  while((up - low) > 1)
    {
    unsigned int mid = (low+up) >> 1;
    if( x > xx[mid]) low = mid;
    else up = mid;
    }
  return low;
}

void interp(SignalType & x, SignalType & xp,
	    SignalType & y, SignalType & result)
{
  result.set_size(x.size());

  unsigned int N = x.size()-1;
  for(unsigned int i = 0; i < result.size(); ++i)
    {
    if( xp[i] <= x[0] ) {result[i] = y[0]; continue;}
    if( xp[i] >= x[N] ) {result[i] = y[N]; continue;}
    unsigned int l = locate(x, xp[i]);
    unsigned int h = l+1;
    double xl = x[l];
    double xh = x[h];
    double yl = y[l];
    double yh = y[h];
    result[i] = yl + (yh-yl)*(xp[i] - xl)/(xh-xl);
    }
}

void print(SignalType & x)
{
  for(unsigned int i = 0; i < x.size(); ++i)
    cout << x[i] << " ";
  cout << endl;
}

void write(const char * fname, SignalType & x)
{
  ofstream ofs(fname);

  for(unsigned int i = 0; i < x.size(); ++i)
    ofs << x[i] << "\n";

  ofs.close();
}

double max(double a, double b)
{
  return (a < b) ? b : a;
}

double csrbf(double r)
{
  double temp = max(0, 1-r);
  return temp*temp*temp*(3*r+1);
}

double phi(double x, vnl_vector<double> &alpha,
	   vnl_vector<double> &beta,
	   vnl_vector<double> &gamma)
{
  double sum = 0;
  for(unsigned int j = 0; j < alpha.size(); ++ j)
    {
    double delta = x - gamma[j];
    sum += alpha[j]*csrbf(sqrt(beta[j]*delta*beta[j]*delta + 0.001));
    }
  return sum;
}

double heaviside(double x)
{
  return 0.5*(1 + (2/PI)*atan(PI*x/EPS));
}

void extract(vnl_vector<double> const &source,
	     unsigned int start,
	     unsigned int end,
	     vnl_vector<double> & result)
{
  unsigned int extractSize = end - start + 1;
  result.set_size(extractSize);
  for(unsigned int i = 0; i < extractSize; ++i)
    {
    result[i] = source[start + i];
    }
}

// result[start:end] = value
void insert(double value,
	    unsigned int start,
	    unsigned int end,
	    vnl_vector<double> & result)
{
  for(unsigned int i = start; i <= end; ++i)
    {
    result[i] = value;
    }
}

double trapz(vnl_vector<double> const & x, double step)
{
  double sum = x[0];
  for(unsigned int i = 1; i < x.size()-1; ++i)
    {
    sum += 2*x[i];
    }
  sum += x[x.size()-1];

  return step*sum/2;
}

class Cost: public vnl_cost_function
{
public:

  Cost(vnl_vector<double> &source,
       vnl_vector<double> &target,
       double weight1,
       double weight2)
    {
      dim = 22;
      I0 = source;
      I1 = target;
      lam = weight1;
      mu = weight2;
    }

  double f(vnl_vector< double > const &theta)
    {
      vnl_vector<double> x;
      double step;
      linspace(x, step, -1, 1, 1000);

      vnl_vector<double> xp = x - theta[0];
      vnl_vector<double> moving;
      interp(x, xp, I0, moving);

      vnl_vector<double> alpha;
      vnl_vector<double> beta;
      vnl_vector<double> gamma;
      extract(theta, 1, 7, alpha);
      extract(theta, 8, 14, beta);
      extract(theta, 15, 21, gamma);

      vnl_vector<double> integrand(x.size());
      for(unsigned int i = 0; i < x.size(); ++i)
      	{
      	double residual = I1[i] - moving[i];
      	integrand[i] = lam*heaviside(-phi(x[i], alpha, beta, gamma))*residual*residual;
	integrand[i] += mu*heaviside(phi(x[i], alpha, beta, gamma));
      	}

      return trapz(integrand, step);
    }

  void gradf(vnl_vector< double > const &theta, vnl_vector< double > &gradient)
    {
      double eps = 1.4901e-08;
      vnl_vector<double> dtheta = theta;

      double current = this->f(theta);
      gradient.set_size(22);
      for(unsigned int i = 0; i < 22; ++i)
	{
	dtheta[i] = theta[i] + eps;
	double next = this->f(dtheta);
	dtheta[i] = theta[i];
	gradient[i] = (next - current)/eps;
	}
    }

  void compute(vnl_vector< double > const &theta, double *f, vnl_vector< double > *g)
    {
      *f = this->f(theta);
      this->gradf(theta, *g);
    }

private:
  Cost(){};

  double lam, mu;
  vnl_vector<double> I0, I1;
};

void kernel(double sigma, vnl_vector<double> &h)
{
  unsigned int sz = 6*ceil(sigma) + 1;
  vnl_vector<double> x;
  double step;
  linspace(x, step, -3*sigma, 3*sigma, sz);
  h.set_size(sz);
  double multiplier = (1/(sqrt(2*PI)*sigma));
  for(unsigned int i = 0; i < sz; ++i)
    {
    h[i] = multiplier*exp(-(x[i]*x[i])/(2*sigma*sigma));
    }
}

void convolve(vnl_vector<double> &signal,
	      vnl_vector<double> &filter,
	      vnl_vector<double> &result)
{
  unsigned int N = signal.size();
  unsigned int M = filter.size();
  unsigned int L = N + M - 1;
  vnl_vector<double> fullresult;
  fullresult.set_size(L);

  // perform convolution
  for(unsigned int l = 0; l < L; ++l)
    {
    double sum = 0;
    for(unsigned int m = 0; m < M; ++m)
      {
      int n = l-m;
      if( (n < 0) || (n >= N) ) continue;
      sum += filter[m]*signal[n];
      }
    fullresult[l] = sum;
    }

  // truncate to same size as input signal
  result.set_size(N);
  unsigned int offset = floor((M-1)/2);
  for(unsigned int n = 0; n < N; ++n)
    {
    result[n] = fullresult[n+offset];
    }
}

void makeI0(vnl_vector<double> & result, unsigned int offset)
{
  vnl_vector<double> square;
  square.set_size(1000);
  square.fill(0);
  insert(1, 250+offset, 750+offset, square);
  vnl_vector<double> h;
  kernel(10,h);
  convolve(square, h, result);
}

void makeI1(vnl_vector<double> & result)
{
  vnl_vector<double> square;
  square.set_size(1000);
  square.fill(0);
  insert(1, 250, 750, square);
  vnl_vector<double> h;
  kernel(10,h);
  convolve(square, h, result);
}

void makeI2(vnl_vector<double> & result)
{
  vnl_vector<double> square;
  square.set_size(1000);
  square.fill(0);
  insert(1, 250, 750, square);
  insert(1.2, 490-50, 510-50, square);
  insert(1.2, 490+200, 510+200, square);
  vnl_vector<double> h;
  kernel(10,h);
  convolve(square, h, result);
}

bool read(SignalType x, string filename)
{
  ifstream infile( filename.c_str() );
  if(!infile) return false;

  bool status = x.read_ascii(infile);
  infile.close();
  return status;
}

int main(int argc, char*argv[])
{
  vul_arg<string> target(0, "Target Signal File");
  vul_arg<std::string> source(0, "Source Signal File");
  vul_arg_parse(argc, argv);

  SignalType T;
  bool ok = read(T, target());
  if(!ok) return EXIT_FAILURE;

  SignalType M;
  makeI0(M, 100);

  vnl_vector<double> initialParameters(22);
  initialParameters[0] = 0; // translation
  initialParameters[1] = -0.2; //initial alpha_1
  initialParameters[2] = 0.2; //initial alpha_2
  initialParameters[3] = -0.2; //initial alpha_3
  initialParameters[4] = 0.2; //initial alpha_4
  initialParameters[5] = -0.2; //initial alpha_5
  initialParameters[6] = 0.2; //initial alpha_6
  initialParameters[7] = -0.2; //initial alpha_7

  insert(2.0, 8, 14, initialParameters); //initial beta_j

  // initial gamma_j
  vnl_vector<double> gamma;
  double step;
  linspace(gamma, step, -1, 1, 7);
  for(unsigned int i = 15, j = 0; i < 23; ++i, ++j) initialParameters[i] = gamma[j];

  print(initialParameters);

  Cost c(M, T, 10, 0.01);

  vnl_lbfgs optimizer(c);

  optimizer.minimize(initialParameters);

  print(initialParameters);
  write("result.dat", initialParameters);

  return 0;
}
