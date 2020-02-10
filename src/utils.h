#pragma once
#include "simplex.h"
#include <dlib/matrix.h>

//To generate cost matrix
matrix<float> gen_K(int m, float eps);
