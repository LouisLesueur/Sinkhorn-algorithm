#pragma once
#include "simplex.h"
#include <dlib/matrix.h>

using namespace std;

//To generate cost matrix
matrix<float> gen_K(int m, float eps);
//The regular Sinkhorn algorithm
simplex bar(const matrix<float> &K, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name);
