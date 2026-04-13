#include "Image.h"


double u_exact(double x, double y, double t, double* D)
{
	double result = 0.0;
	double det = D[0] * D[3] - D[1] * D[2];
	std::vector<std::vector<double>> A = {
		{D[3] / det, -D[1] / det},
		{-D[2] / det, D[0] / det}
	};

	double koeff = 1. / (4. * M_PI * det * t);
	double exp_arg = -(A[0][0] * x * x + (A[0][1] + A[1][0]) * x * y + A[1][1] * y * y) / (4 * t);
	result = koeff * exp(exp_arg);
	return result;
}

Image::Image(uchar* data, int w, int h, int bytesPerLine)
{
	width = w;
	height = h;
	size = w * h;

	imageData.resize(size);
	damaged.resize(size);
	mask.resize(size);
	laplace.resize(size);
	smoothed.resize(size);

	for (int i = 0; i < height; i++) {
		const uchar* row = data + i * bytesPerLine;
		int base = i * width;
		for (int j = 0; j < width; j++) {
			imageData[base + j] = row[j] / 255.0;
		}
	}
}

void Image::generateMask(int p)
{
	int toRemove = static_cast<int>(size * p / 100.0 + 0.5);

	std::fill(mask.begin(), mask.end(), 1);

	std::vector<int> indexy(size);
	std::iota(indexy.begin(), indexy.end(), 0);		//toto naplna indexy postupne 0,1,2,...,size-1

	static std::mt19937 rng{ std::random_device{}() };	//random generator
	std::shuffle(indexy.begin(), indexy.end(), rng); 	//shuffle nam zamiesa indexy, aby sme mohli z nich vybrat nahodne

	for (int i = 0; i < toRemove; i++) {
		mask[indexy[i]] = 0;
	}

	for (int i = 0; i < size; i++) {
		damaged[i] = mask[i] ? imageData[i] : 0.0;
	}

	int removedCount = 0;
	for (int i = 0; i < size; i++) {
		if (mask[i] == 0)
			removedCount++;
	}

	std::cout << "Removed pixels: "
		<< removedCount << " / " << size
		<< " (" << (100.0 * removedCount / size) << "%)\n";

}

void Image::Laplace()
{
	//konstrukcia Laplaceovej matice

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size); Eigen::VectorXd b = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd xs(size);
	M.reserve(Eigen::VectorXi::Constant(size, 5));

	for (int i = 0; i < height; i++) 
	{
		int row = i * width;
		for (int j = 0; j < width; j++) 
		{
			int id = row + j;
			//odstraneny pixel
			if (mask[id] == 0) 
			{
				//ROHY
				if (i == 0 && j == 0)	//vlavo dole
				{
					M.insert(id, id) = 4.;
					M.insert(id, id + width) = -2.0;
					M.insert(id, id + 1) = -2.0;
				}
				else if (i == 0 && j == width - 1)	//vpravo dole
				{
					M.insert(id, id) = 4.;
					M.insert(id, id + width) = -2.0;
					M.insert(id, id - 1) = -2.0;
				}
				else if (i == height - 1 && j == 0)	//vlavo hore
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -2.0;
					M.insert(id, id + 1) = -2.0;
				}
				else if (i == height - 1 && j == width - 1)	//vpravo hore
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -2.0;
					M.insert(id, id - 1) = -2.0;
				}

				//HRANY
				else if (i == 0 && j>=1 && j<=width-2)	//dole
				{
					M.insert(id, id) = 4.;
					M.insert(id, id + width) = -2.0;
					M.insert(id, id - 1) = -1.0;
					M.insert(id, id + 1) = -1.0;
				}
				else if (i == height - 1 && j>=1&&j<=width-2) //hore
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -2.0;
					M.insert(id, id - 1) = -1.0;
					M.insert(id, id + 1) = -1.0;
				}
				else if (j == 0 && i>=1 && i<=height-2)	//vlavo
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -1.0;
					M.insert(id, id + width) = -1.0;
					M.insert(id, id + 1) = -2.0;
				}
				else if (j == width - 1 && i >= 1 && i <= height - 2)	//vpravo
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -1.0;
					M.insert(id, id + width) = -1.0;
					M.insert(id, id - 1) = -2.0;
				}
				//ZVYSOK
				else
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -1.0;
					M.insert(id, id + width) = -1.0;
					M.insert(id, id - 1) = -1.0;
					M.insert(id, id + 1) = -1.0;
				}
			}
			else 
			{
				M.insert(id, id) = 1.0;
				b(id) = damaged[id];
			}
		}
	}

	//SOLVE
	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

	solver.analyzePattern(M); solver.factorize(M);
	if (solver.info() != Eigen::Success) {
		std::cout << "Error in factorization of matrix" << std::endl;
		return;
	}
	xs = solver.solve(b);
	if (solver.info() != Eigen::Success) 
	{
		std::cout << "Error in solver" << std::endl;
	}
	for (int i = 0; i < size; i++) 
	{
		laplace[i] = xs[i];
	}


	
}

void Image::Smooth(double lambda)
{
	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size); Eigen::VectorXd b = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd xs(size);
	M.reserve(Eigen::VectorXi::Constant(size, 5));

	for (int i = 0; i < height; i++)
	{
		int row = i * width;
		for (int j = 0; j < width; j++)
		{
			int id = row + j;

			if (i == 0 && j == 0)	//vlavo dole
			{
				M.insert(id, id) = 4. +lambda;
				M.insert(id, id + width) = -2.0;
				M.insert(id, id + 1) = -2.0;
			}
			else if (i == 0 && j == width - 1)	//vpravo dole
			{
				M.insert(id, id) = 4. +lambda;
				M.insert(id, id + width) = -2.0;
				M.insert(id, id - 1) = -2.0;
			}
			else if (i == height - 1 && j == 0)	//vlavo hore
			{
				M.insert(id, id) = 4.+lambda;
				M.insert(id, id - width) = -2.0;
				M.insert(id, id + 1) = -2.0;
			}
			else if (i == height - 1 && j == width - 1)	//vpravo hore
			{
				M.insert(id, id) = 4.+lambda;
				M.insert(id, id - width) = -2.0;
				M.insert(id, id - 1) = -2.0;
			}

			//HRANY
			else if (i == 0 && j >= 1 && j <= width - 2)	//dole
			{
				M.insert(id, id) = 4. +lambda;
				M.insert(id, id + width) = -2.0;
				M.insert(id, id - 1) = -1.0;
				M.insert(id, id + 1) = -1.0;
			}
			else if (i == height - 1 && j >= 1 && j <= width - 2) //hore
			{
				M.insert(id, id) = 4.+lambda;
				M.insert(id, id - width) = -2.0;
				M.insert(id, id - 1) = -1.0;
				M.insert(id, id + 1) = -1.0;
			}
			else if (j == 0 && i >= 1 && i <= height - 2)	//vlavo
			{
				M.insert(id, id) = 4. +lambda;
				M.insert(id, id - width) = -1.0;
				M.insert(id, id + width) = -1.0;
				M.insert(id, id + 1) = -2.0;
			}
			else if (j == width - 1 && i >= 1 && i <= height - 2)	//vpravo
			{
				M.insert(id, id) = 4. +lambda;
				M.insert(id, id - width) = -1.0;
				M.insert(id, id + width) = -1.0;
				M.insert(id, id - 1) = -2.0;
			}
			//ZVYSOK
			else
			{
				M.insert(id, id) = 4. +lambda;
				M.insert(id, id - width) = -1.0;
				M.insert(id, id + width) = -1.0;
				M.insert(id, id - 1) = -1.0;
				M.insert(id, id + 1) = -1.0;
			}
			b(id) = lambda*laplace[id];
		}
	}

	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

	solver.analyzePattern(M); solver.factorize(M);
	if (solver.info() != Eigen::Success) {
		std::cout << "Error in factorization of matrix" << std::endl;
		return;
	}
	xs = solver.solve(b);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "Error in solver" << std::endl;
	}
	for (int i = 0; i < size; i++)
	{
		smoothed[i] = xs[i];
	}


}
double Image::Anisotropic(bool exportData)
{
	setD();

	int width = N + 2;
	int height = N + 2;
	int size = width * height;

	std::vector<double> u_old(size, 0.0);

	double h = 2.0 / N;
	double tau = h * h;
	double Tstart = 0.2;
	double Tend = 0.3;
	double error_sum = 0.0;

	// ---------------- INIT ----------------
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}

	// initial error
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 + j * h;
			int k = i * width + j;

			double exact = u_exact(x, y, Tstart, &D[0][0]);
			double diff = u_old[k] - exact;
			error_sum += diff * diff;
		}
	}

	// ---------------- EXPORT INIT ----------------
	std::string folder;
	if (exportData)
	{
		folder = "anisotropic_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta));

		std::filesystem::create_directories(folder);

		std::ofstream file(folder + "/step_0.csv");
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				file << u_old[i * width + j];
				if (j < width - 1) file << ",";
			}
			file << "\n";
		}
		file.close();
	}

	// ---------------- BUILD MATRIX ONCE ----------------
	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;

			bool isBoundary = (i == 0 || i == height - 1 || j == 0 || j == width - 1);

			if (isBoundary)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			int kE = k + 1;
			int kW = k - 1;
			int kN = k + width;
			int kS = k - width;

			int kNE = kN + 1;
			int kNW = kN - 1;
			int kSE = kS + 1;
			int kSW = kS - 1;

			double bE = alpha;
			double bW = alpha;
			double bN = beta;
			double bS = beta;
			double bNE = gamma / 2.0;
			double bNW = -gamma / 2.0;
			double bSW = gamma / 2.0;
			double bSE = -gamma / 2.0;

			double bcenter = (bE + bW + bN + bS + bNE + bNW + bSE + bSW);

			M.insert(k, k) = 1.0 + scale * bcenter;

			M.insert(k, kE) = -scale * bE;
			M.insert(k, kW) = -scale * bW;
			M.insert(k, kN) = -scale * bN;
			M.insert(k, kS) = -scale * bS;

			M.insert(k, kNE) = -scale * bNE;
			M.insert(k, kNW) = -scale * bNW;
			M.insert(k, kSE) = -scale * bSE;
			M.insert(k, kSW) = -scale * bSW;
		}
	}

	M.makeCompressed();

	// ---------------- FACTORIZE ONCE ----------------
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);

	// ---------------- TIME LOOP ----------------
	int steps = static_cast<int>((Tend - Tstart) / tau);

	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		Eigen::VectorXd b(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 + j * h;

				bool isBoundary = (i == 0 || i == height - 1 || j == 0 || j == width - 1);

				if (isBoundary)
				{
					b(k) = u_exact(x, y, t_n + tau, &D[0][0]);
				}
				else
				{
					b(k) = u_old[k];
				}
			}
		}

		Eigen::VectorXd sol = solver.solve(b);

		for (int k = 0; k < size; k++)
			u_old[k] = sol(k);

		// accumulate error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0 + j * h;
				int k = i * width + j;

				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;

				error_sum += diff * diff;
			}
		}
	}

	// ---------------- EXPORT FINAL ----------------
	if (exportData)
	{
		std::ofstream finalFile(folder + "/final.csv");

		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				finalFile << u_old[i * width + j];
				if (j < width - 1) finalFile << ",";
			}
			finalFile << "\n";
		}
		finalFile.close();
	}

	return std::sqrt(tau * h * h * error_sum);
}
void Image::setD()
{
	if (theta == 10.0) {
		D[0][0] = 0.488;
		D[0][1] = -0.086;
		D[1][0] = -0.086;
		D[1][1] = 0.016;
	}
	else if (theta == 22.5) {
		D[0][0] = 0.431;
		D[0][1] = -0.178;
		D[1][0] = -0.178;
		D[1][1] = 0.074;
	}
	else if (theta == 45.0) {
		D[0][0] = 0.253;
		D[0][1] = -0.252;
		D[1][0] = -0.252;
		D[1][1] = 0.253;
	}
}

