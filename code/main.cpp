#include "simplex.h"
#include "pi.h"

using namespace std;

int main()
{
	InitRandom();
	simplex s1(100);
	simplex s2(100);
	double eps = 0.01;
	int n_iter = 10;
	pi gamma = W(s1, s2, eps, n_iter);
	gamma.plot();

	return 0;
}
