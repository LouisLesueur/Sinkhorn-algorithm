#pragma once
#include "simplex.h"
#include <dlib/matrix.h>

using namespace std;

matrix<float> gen_K(int m, float eps);
matrix<float> gen_C(int m);
simplex bar(const matrix<float> &K, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name);
simplex bar_log(const matrix<float> &C, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name);
