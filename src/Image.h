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

	//toto bude na ukladanie time stepov toho posledneho zadania
	std::vector<std::vector<double>> evolutionField1;
	std::vector<std::vector<double>> evolutionField2;
	std::vector<std::vector<double>> evolutionField3; 


	//ANIZOTROPNE
	int N = 0;
	double theta = 0.0;
	double D[2][2];

	//posledne zadanie
	int timeSteps = 0;
	double tau = 0.0;
	double K1 = 0.0;
	double K2 = 0.0;


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
	double* getEvolutionFrame(int fieldId, int step);
	int getEvolutionFrameCount(int fieldId) const;

	void generateMask(int p);
	void Laplace();
	void Smooth(double lambda);

	void setN (int n) { N = n; }
	void setTheta(double t) { theta = t; }
	double Anisotropic_Classic(bool exportData);	//klasicka anizotropna difuzia
	double Anisotropic_Modified(bool exportData);	//modifikovana anizotropna difuzia (ADCM)
	double S1_FBDS_Classic(bool exportData);		//S1 Forward-Backward Diffusion Scheme - klasicka anizotropna difuzia
	double S2_FBDS_Classic(bool exportData);		//S2 Forward-Backward Diffusion Scheme - klasicka anizotropna difuzia
	double S1_FBDS_ADCM(bool exportData);			//S1 Forward-Backward Diffusion Scheme - modifikovana anizotropna difuzia (ADCM)
	double S2_FBDS_ADCM(bool exportData);			//S2 Forward-Backward Diffusion Scheme - modifikovana anizotropna difuzia (ADCM)


	void setD();

	void generateRandomImage(int n);
	void setTimeSteps(int steps) { timeSteps = steps; }
	void setTau(double t) { tau = t; }
	void setK1(double k) { K1 = k; }
	void setK2(double k) { K2 = k; }
	void variableDCM_forField(int fieldId, std::vector<std::vector<double>>& storage);
	void variableDCM();
	void computeD(double v1, double v2, double D[2][2]);	//funkcia na vypocet tych specifickych D podla vektoroveho pola
	void vectorField(double x1, double x2, double& v1, double& v2, int fieldId);	//tu sa bude pocitat to vektorove pole samotne
};