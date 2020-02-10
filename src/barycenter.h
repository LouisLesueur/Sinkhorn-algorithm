#pragma once
#include "simplex.h"
#include "utils.h"
#include <dlib/matrix.h>

using namespace std;

//The regular Sinkhorn algorithm
simplex bar(const matrix<float> &K, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name);
