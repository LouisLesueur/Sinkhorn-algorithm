#pragma once
#include "simplex.h"
#include<vector>

using namespace std;

simplex bar(const simplex & p1, const simplex & p2, double lambda, double eps, int n_iter, string name);