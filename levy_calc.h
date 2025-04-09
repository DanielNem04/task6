#include <cmath>
#include <complex>

using std::complex;

#define NMIN_SIN 60 // at least this much points in a quarter period of the sine function
#define NMIN_EXP 50 // at least this much points in one 1/e change of the exponential
#define LIMIT 30    // if the exponent in the exponential in the iterand is less than minus LIMIT, stop iteration

#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// for the Lanczos approximation of the gamma function
const double lanczos_coeff[9] = {
0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
const double Euler_gamma = 0.57721566490153286060;
const double g_lanczos = 7.;
const complex<double> I(0., 1.);

complex<double> Gamma(complex<double> z)  //Utilizing the "Lanczos approximation"
{
  complex<double> A_g = lanczos_coeff[0]; 
  for(int i=1; i <= (int)(g_lanczos + 1.); i++)
    A_g += lanczos_coeff[i] / (z - 1. + double(i));
  return sqrt(2. * M_PI) * pow( z + g_lanczos - 0.5, z - 0.5 ) * exp(- (z + g_lanczos - 0.5)) * A_g;
}

double Gamma(double x)
{
  complex<double> z(x, 0.);
  return real(Gamma(z));
}

double levy_calc(double x, double R, double alpha)
{
  double result = 0.;  
  double t = 0.;
  if(x < 0. || alpha < 0.5) return 0;
  if(x < 0.001)
    return  pow(2., 3 / alpha - 1.) / alpha / SQR(M_PI) / CUBE(R) * Gamma(3. / alpha);
  
  double dt_sin = M_PI / 2. / NMIN_SIN;
  double dt_exp = pow(2., 1.0 / alpha) * x / NMIN_EXP;
  double dt = (dt_sin < dt_exp) ? dt_sin : dt_exp;
  double exponent_dt = 0., exponent_mean = 0.;
  double exponential = 1., exponential_dt = 1., exponential_mean = 1.;
  double sinust = 0., sinust_dt = 0., sinust_mean = 0.;
  for(;;)
  {
    exponent_dt = -pow(t + dt, alpha ) / ( 2 * pow(x, alpha) );
    exponent_mean =  -pow((2 * t + dt) / 2, alpha ) / ( 2 * pow(x, alpha) );
    exponential = exponential_dt;
    exponential_mean = exp(exponent_mean);
    exponential_dt = exp(exponent_dt);
    sinust = sinust_dt;
    sinust_mean = sin((2. * t + dt) / 2.);
    sinust_dt = sin(t + dt);
    if(exponent_dt + log(1. + t) < -LIMIT) // change was here
      break;
    double parabolic_sum = (dt / 6.) * ( (exponential * sinust * t) + 4. * (exponential_mean * sinust_mean * ((2. * t + dt) / 2.)) + (exponential_dt * sinust_dt * (t + dt)));
    result += parabolic_sum;
    t += dt; 
  }
  return result * 1.0/(2. * SQR(M_PI) * CUBE(R * x));
}
