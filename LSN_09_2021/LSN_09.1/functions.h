#include "TSP_position.h"
#include <vector>
#include <algorithm>

using namespace std;

double error(std::vector <double> v1, std::vector  <double> v2, int n);
void data_blocking(std::vector  <double> value, std::vector  <double> &sum_prog, std::vector  <double> &err_prog);
double TSP_cost(vector <int> o, TSP_position p[]);
bool TSP_check(std::vector <int> o);

void crossover(vector <int> &v1, vector <int> &v2, int cut);
double expo(double y, double lambda);
int selection(vector <double> fitness, double rand);

//Mutations
void swap (std::vector <int> &vett, int a, int b);
void shift (vector <int> &vett, int shift);
void invert_order(vector <int> &vett, int a, int b);
void part_perm(vector <int> &vett, int a, int b, int c, int d);


