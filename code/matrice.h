#pragma once
#include <cassert>
#include <iostream>

using namespace std;


class Matrice {

    private:
        //m lig et n col
        int m,n;
	double* tab;
        int* compteur;

    public:
        int nlignes() const {return m;}
        int ncolonnes() const {return n;}

        //Opérateurs
        double operator()(int i, int j) const;
        double& operator()(int i, int j);

	Matrice operator*(double k);
	Matrice operator/(double k);

	Matrice &operator=(const Matrice &);

        //Constructeurs/destructeur
        Matrice(int nligcol);
        Matrice(int nlig, int ncol);
        Matrice(const Matrice &A);
        ~Matrice();

};

Matrice operator*(const Matrice& A, const Matrice& B);
Matrice operator/(const Matrice& A, const Matrice& B);
ostream& operator<<( ostream& str, const Matrice& A);
Matrice transpose(Matrice A);
