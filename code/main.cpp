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
	cout << "Ok ça c'était cool mais passons aux choses sérieuses" << endl;
	simplex s1(4);
	simplex s2(4);
	double eps = 0.01;
	int n_iter = 10;
	pi gamma = W(s1, s2, eps, n_iter);
	cout << "Premier simplex " << s1 << endl;
	cout << "Premiere marginale " << gamma.first_marginal() << endl;
	cout << "Second simplex " << s2 << endl;
	cout << "Seconde marginale " << gamma.second_marginal() << endl;

	return 0;
}
