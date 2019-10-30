#include "simplex.h"
#include "pi.h"

using namespace std;

int main()
{
	InitRandom();
	simplex b(4);
	std::cout<<"test simplexe: "<<b<<std::endl;
	pi p(4);
	std::cout<<"test pi: "<<std::endl;
	std::cout<<"--première marginale: "<<p.first_marginal()<<std::endl;
	std::cout<<"--seconde  marginale: "<<p.second_marginal()<<std::endl;
	return 0;
}
