#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <dlib/gui_widgets.h>
#include <dlib/pixel.h>
#include <dlib/matrix.h>
#include <dlib/array2d.h>
#include <dlib/image_io.h>

using namespace dlib;
using namespace std;

class simplex{
	private:
		matrix<float> tab;
		string name;
		int constant, width, height;
	public:
		//empty constructor
		simplex(){}
		//build a simplex from existing data
		simplex(matrix<float> VAL, int WIDTH, int HEIGHT, int CSTE, string path);
		//To build a simplex from a image
		simplex(string path);

		//To get the characteritics of the simplex
		int length() const{return tab.nr()*tab.nc();};
		int w() const{return width;};
		int h() const{return height;};
		int cte() const{return constant;};
		matrix<float> val() const{return tab;};

		//To export as an image
		void export_to_img();

		//To gets elements of the simplex
		float operator()(int i) const;
		float& operator()(int i);

};
