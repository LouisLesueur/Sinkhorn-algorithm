#include <Eigen/Dense>
#include <Eigen/Sparse>
#undef Success  
#include "pi.h"
#include <vector>
#include <cmath>

using namespace Eigen;
using namespace std;

//===================================CLASSE PI=========================

//-------------------------CONSTRUCTEUR--------------------------------
pi::pi(int n){
	size = n;
	tab = new double[n*n];
	double sum = 0;

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			double coef = Random(0,100);
			tab[i+n*j] = coef;
			sum += coef;
		}
	}
	for(int i=0; i<n*n; i++)
		tab[i] /= sum;
}

//----------------------------OPERATEURS----------------------------------

double pi::operator()(int i, int j) const{
	assert(i<size && j<size);
	return tab[i+size*j];
}

double& pi::operator()(int i, int j){
	assert(i<size && j<size);
	return tab[i+size*j];
}

//------------------------MARGINALES---------------------------------------

simplex pi::first_marginal(){
	double marg[size];
	for(int i=0; i<size; i++){
		double sum = 0;
		for(int j=0; j<size; j++)
			sum += tab[i+size*j];
		marg[i] = sum;
	}
	return simplex(marg, size, '1');
}

simplex pi::second_marginal(){
	double marg[size];
	for(int j=0; j<size; j++){
		double sum = 0;
		for(int i=0; i<size; i++)
			sum += tab[i+size*j];
		marg[j] = sum;
	}
	return simplex(marg, size, '2');
}

//====================================================================

void pi::plot(){
	ofstream monFlux("pi.csv");
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++)
			monFlux<<" "<<double(tab[i+size*j]);
		monFlux<<endl;
	}
}


SparseMatrix<double> div(const SparseMatrix<double> a, const SparseMatrix<double> b){
	int n = a.rows();
	int m = a.cols();
	SparseMatrix<double> c(n,m);
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			if(b.coeff(i,j) == 0)
				cout << "Division par 0 !" << endl;
		c.insert(i,j) = a.coeff(i,j)/b.coeff(i,j);
		}
	}
	return c;
}

SparseMatrix<double> product(const SparseMatrix<double>& A, const SparseMatrix<double>& B){
	SparseMatrix<double> C(A.rows(), A.cols());
		for(int i=0; i<A.rows(); i++){
			for(int j=0; j<A.cols(); j++){
				C.insert(i,j) = A.coeff(i,j)*B.coeff(i,j);
			}
		}
	return C;
}

SparseMatrix<double> matrix_pow(const SparseMatrix<double>& A, double lambda){
	SparseMatrix<double> C(A.rows(), A.cols());
	for(int i=0; i<A.rows(); i++){
		for(int j=0; j<A.cols(); j++){
			C.insert(i,j) = pow(A.coeff(i,j), lambda);
		}
	}
	return C;
}



pi W(const simplex &s1, const simplex &s2, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	SparseMatrix<double> ksi(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			ksi.insert(i,j) = exp(-dis/eps);
		}
	}
	SparseVector<double> p1(m);
	SparseVector<double> p2(m);
	SparseVector<double> v(m);
	for(int i=0; i<m; i++){
		p1.insert(i) = s1(i);
		p2.insert(i) = s2(i);
		v.insert(i) = 1;
	}
	// Using recursion formula,
	for(int i=0; i<n_iter; i++){
		v = div(p2, ksi.transpose()*div(p1, ksi*v));
	}

	// Creating diag matrix
	SparseVector<double> u(m);
	u = div(p1, ksi*v);
	SparseMatrix<double> diag1(m,m);
	SparseMatrix<double> diag2(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			if(i==j){
				diag1.insert(i,j) = u.coeff(i);
				diag2.insert(i,j) = v.coeff(i);
			}
			else{
				diag1.insert(i,j) = 0;
				diag2.insert(i,j) = 0;
			}
		}
	}
	SparseMatrix<double> final_prod = diag1*ksi*diag2;
	pi gamma(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			gamma(i,j) = final_prod.coeff(i,j);
		}
	}
	return gamma;
}

simplex bar(const simplex & s1, const simplex & s2, double lambda, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	SparseMatrix<double> ksi(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			ksi.insert(i,j) = exp(-dis/eps);
		}
	}
	SparseVector<double> p1(m), p2(m),
		 v1(m), v2(m),
		 u1(m), u2(m);

	for(int i=0; i<m; i++){ // ukn and vkn
		p1.insert(i) = s1(i);
		p2.insert(i) = s2(i);
		v1.insert(i) = 1;
		v2.insert(i) = 1;
		u1.insert(i) = 1;
		u2.insert(i) = 1;
	}
	for(int i=0; i<n_iter; i++){
		// Computing vk(n+1)
		v1 = div(p1, ksi.transpose()*u1);
		v2 = div(p1, ksi.transpose()*u2);
		// Computing p(n+1)
		SparseMatrix<double> p = product(matrix_pow((product(u1, ksi*v1)), lambda), matrix_pow(product(u2, ksi*v2), 1-lambda));
		// Computing uk(n+1)
		u1 = div(p, ksi*v1);
		u2 = div(p, ksi*v2);
	}
	// Computing the first projection matrix
	SparseMatrix<double> diag1(m,m);
	SparseMatrix<double> diag2(m,m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			if(i==j){
				diag1.insert(i,j) = u1.coeff(i);
				diag2.insert(i,j) = v1.coeff(i);
			}
			else{
				diag1.insert(i,j) = 0;
				diag2.insert(i,j) = 0;
			}
		}
	}
	SparseMatrix<double> final_prod = diag1*ksi*diag2;
	pi gamma(m);
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			gamma(i,j) = final_prod.coeff(i,j);
	simplex first_marginal = gamma.first_marginal();
	simplex barycentre(m, '3');
	for(int i=0; i<m; i++)
		barycentre(i) = first_marginal(i);
	return barycentre;
}
