#pragma once
#include <cassert>
#include <iostream>
#include <ctime>

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

        //Constructeurs/destructeur
        Matrice(int nligcol);
        Matrice(int nlig, int ncol);
        Matrice(const Matrice &A);
        ~Matrice();

};

Matrice operator*(const Matrice& A, const Matrice& B);
ostream& operator<<( ostream& str, const Matrice& A);
