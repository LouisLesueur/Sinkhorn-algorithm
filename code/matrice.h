#pragma once

class Matrice{
public:
	void Affiche(); 
	double get(int i, int j) const; // Accesseur
	void set(double x, int i, int j); // Mutateur
	int nlig; // Nb lignes
	int ncol; // Nb colonnes
	Matrice operator*(Matrice & M);
	double operator()(int i, int j); // "Accesseur amélioré"
	void Alea(); // Attribue des coefficients aléatoires à la matrice	
	Matrice(int n, int m); // Constructeur
	Matrice(const Matrice & A); // Constructeur par copie
	~Matrice(); // Destructeur
private:
	double* tab;
	int* compteur;
};

