#pragma once
#include <QVector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

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

public:
	Image() {};
	Image(uchar* data, int w, int h);

	int getwidth() { return width; }
	int getheight() { return height; }

	double* getImageData() { return imageData.data(); }
	uint8_t* getMask() { return mask.data(); }
	double* getDamaged() { return damaged.data(); }
	double* getLaplace() { return laplace.data(); }

	void generateMask(int p);
};