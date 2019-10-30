#pragma once

#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>

//À appeler une fois dans le main
void InitRandom();
//Pour générer un double aléatoire entre a et b
double Random(double a, double b);

class simplex{
	private:
		double* tab;
		int size;
	public:
		simplex(int n);
		int length() const{return size;};

		double operator()(int i) const;
		double& operator()(int i);

		~simplex(){delete [] tab;}

};

std::ostream& operator<<(std::ostream& str, const simplex& s);
