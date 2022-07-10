#include <vector>

double error(std::vector <double> v1, std::vector  <double> v2, int n);
void data_blocking(std::vector  <double> value, std::vector  <double> &sum_prog, std::vector  <double> &err_prog);
double psi2_100(double x, double y, double z);
double psi2_210(double x, double y, double z);
//void Metropolis_uniform(std::vector <double> x, std::vector <double> y, std::vector <double> z, std::vector <double> r, int i, double q_xy, double delta);
//void Metropolis_gauss(std::vector <double> x, std::vector <double> y, std::vector <double> z, std::vector <double> r, int i, double q_xy, double delta);
