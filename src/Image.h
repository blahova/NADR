#pragma once
#define _USE_MATH_DEFINES
#include <QVector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <random>
#include <math.h>

//eigen
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Eigen/SparseLU"

class Image {
private:
	int width = 0;
	int height = 0;
	int size = 0;

	std::vector<double> imageData;	//originalny ovrazok
	std::vector<uint8_t> mask;	//maska
	std::vector<double> damaged;	//aplikacia masky na obrazok
	std::vector<double> laplace;	//uz po laplaceovi
	std::vector<double> smoothed;	//smooth

	//ANIZOTROPNE
	int N = 0;
	double theta = 0.0;
	double D[2][2];


public:
	Image() {};
	Image(uchar* data, int w, int h, int bytesPerLine);

	int getwidth() { return width; }
	int getheight() { return height; }

	double* getImageData() { return imageData.data(); }
	uint8_t* getMask() { return mask.data(); }
	double* getDamaged() { return damaged.data(); }
	double* getLaplace() { return laplace.data(); }
	double* getSmoothed() { return smoothed.data(); }

	void generateMask(int p);
	void Laplace();
	void Smooth(double lambda);

	void setN (int n) { N = n; }
	void setTheta(double t) { theta = t; }
	double Anisotropic(bool exportData);
	void setD();
};