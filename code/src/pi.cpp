#include "pi.h"
#include <vector>
#include <cmath>
#include <dlib/matrix.h>

using namespace dlib;
using namespace std;

//===================================CLASSE PI=========================

//-------------------------CONSTRUCTEUR--------------------------------
PI::PI(int n){
	size = n;
	tab = new double[n*n];
	/*double sum = 0;

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			double coef = Random(0,100);
			tab[i+n*j] = coef;
			sum += coef;
		}
	}
	for(int i=0; i<n*n; i++)
		tab[i] /= sum;*/
}

//----------------------------OPERATEURS----------------------------------

double PI::operator()(int i, int j) const{
	assert(i<size && j<size);
	return tab[i+size*j];
}

double& PI::operator()(int i, int j){
	assert(i<size && j<size);
	return tab[i+size*j];
}

//------------------------MARGINALES---------------------------------------

simplex PI::first_marginal(){
	double marg[size];
	for(int i=0; i<size; i++){
		double sum = 0;
		for(int j=0; j<size; j++)
			sum += tab[i+size*j];
		marg[i] = sum;
	}
	return simplex(marg, size, size, 1, '1');
}

simplex PI::second_marginal(){
	double marg[size];
	for(int j=0; j<size; j++){
		double sum = 0;
		for(int i=0; i<size; i++)
			sum += tab[i+size*j];
		marg[j] = sum;
	}
	return simplex(marg, size, size, 1, '2');
}

//====================================================================

void PI::plot(){
	ofstream monFlux("pi.csv");
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++)
			monFlux<<" "<<double(tab[i+size*j]);
		monFlux<<endl;
	}
}


matrix<double> div(const matrix<double> a, const matrix<double> b){
	int n = a.nr();
	int m = a.nc();
	matrix<double> c(n,m);
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			if(b(i,j) == 0)
				cout << "Division par 0 !" << endl;
		c(i,j) = a(i,j)/b(i,j);
		}
	}
	return c;
}

matrix<double> product(const matrix<double>& A, const matrix<double>& B){
	matrix<double> C(A.nr(), A.nc());
		for(int i=0; i<A.nr(); i++){
			for(int j=0; j<A.nc(); j++){
				C(i,j) = A(i,j)*B(i,j);
			}
		}
	return C;
}

matrix<double> matrix_pow(const matrix<double>& A, double lambda){
	matrix<double> C(A.nr(), A.nc());
	for(int i=0; i<A.nr(); i++){
		for(int j=0; j<A.nc(); j++){
			C(i,j) = pow(A(i,j), lambda);
		}
	}
	return C;
}



PI W(const simplex &s1, const simplex &s2, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	matrix<double> ksi(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			ksi(i,j) = exp(-dis/eps);
		}
	}
	matrix<double> p1(m,1);
	matrix<double> p2(m,1);
	matrix<double> v(m,1);
	for(int i=0; i<m; i++){
		p1(i,0) = s1(i);
		p2(i,0) = s2(i);
		v(i,0) = 1;
	}
	// Using recursion formula,
	for(int i=0; i<n_iter; i++){
		v = div(p2, trans(ksi)*div(p1, ksi*v));
	}

	// Creating diag matrix
	matrix<double> u(m,1);
	u = div(p1, ksi*v);
	matrix<double> diag1(m,m);
	matrix<double> diag2(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			if(i==j){
				diag1(i,j) = u(i,0);
				diag2(i,j) = v(i,0);
			}
			else{
				diag1(i,j) = 0;
				diag2(i,j) = 0;
			}
		}
	}
	matrix<double> final_prod = diag1*ksi*diag2;
	PI gamma(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			gamma(i,j) = final_prod(i,j);
		}
	}
	return gamma;
}

simplex bar(const simplex & s1, const simplex & s2, double lambda, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	matrix<double> ksi(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			ksi(i,j) = exp(-dis/eps);
		}
	}
	matrix<double> p1(m,1), p2(m,1),
		 v1(m,1), v2(m,1),
		 u1(m,1), u2(m,1);

	for(int i=0; i<m; i++){ // ukn and vkn
		p1(i,0) = s1(i);
		p2(i,0) = s2(i);
		v1(i,0) = 1;
		v2(i,0) = 1;
		u1(i,0) = 1;
		u2(i,0) = 1;
	}
	for(int i=0; i<n_iter; i++){
		// Computing vk(n+1)
		v1 = div(p1, trans(ksi)*u1);
		v2 = div(p1, trans(ksi)*u2);
		// Computing p(n+1)
		matrix<double> p = product(matrix_pow((product(u1, ksi*v1)), lambda), matrix_pow(product(u2, ksi*v2), 1-lambda));
		// Computing uk(n+1)
		u1 = div(p, ksi*v1);
		u2 = div(p, ksi*v2);
	}
	// Computing the first projection matrix
	matrix<double> diag1(m,m);
	matrix<double> diag2(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			if(i==j){
				diag1(i,j) = u1(i,0);
				diag2(i,j) = v1(i,0);
			}
			else{
				diag1(i,j) = 0;
				diag2(i,j) = 0;
			}
		}
	}
	matrix<double> final_prod = diag1*ksi*diag2;
	PI gamma(m);
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			gamma(i,j) = final_prod(i,j);
	simplex first_marginal = gamma.first_marginal();
	simplex barycentre(m, '3');
	for(int i=0; i<m; i++)
		barycentre(i) = first_marginal(i);
	return barycentre;
}
