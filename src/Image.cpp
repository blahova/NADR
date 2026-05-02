#include "Image.h"


double u_exact(double x, double y, double t, double* D)
{
	double result = 0.0;
	double det = D[0] * D[3] - D[1] * D[2];
	std::vector<std::vector<double>> A = {
		{D[3] / det, -D[1] / det},
		{-D[2] / det, D[0] / det}
	};

	double koeff = 1. / (4. * M_PI * sqrt(det) * t);
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


double Image::Anisotropic_Classic(bool exportData)
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

	//INIT 
	for (int i = 0; i < height; i++)
	{
		double y = -1.0-h/2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0-h/2 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}

	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	const int d[8] = { +1, -1, +width, -width,
					+width + 1, +width - 1,
					-width + 1, -width - 1 };

	double b[8] = { alpha, alpha,
					beta,  beta,
					gamma / 2.0, -gamma / 2.0,
					-gamma / 2.0, gamma / 2.0 };  // E,W,N,S,NE,NW,SE,SW

	std::string folder = "DCM_Classic";
	if (exportData)
	{
		
		std::filesystem::create_directories(folder);
		std::string filename = folder + "/step_0_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream file(filename);

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

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;
			if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			double bcenter = 0.0;
			for (int s = 0; s < 8; s++)
			{
				int q = k + d[s];
				M.insert(k, q) = -scale * b[s];
				bcenter += b[s];
			}
			M.insert(k, k) = 1.0 + scale * bcenter;
		}
	}

	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);

	int steps = static_cast<int>(std::round((Tend - Tstart) / tau));


	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		Eigen::VectorXd b(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 -h/2. + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 -h/2. + j * h;

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

		// error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 -h/2 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0-h/2 + j * h;
				int k = i * width + j;

				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;

				error_sum += diff * diff;
			}
		}
	}

	// EXPORT FINAL
	if (exportData)
	{
		std::string finalName = folder + "/final_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream finalFile(finalName);

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

double Image::Anisotropic_Modified(bool exportData)
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

	//INIT 
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 - h / 2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 - h / 2 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}

	std::string folder = "ADCM_Modified";
	if (exportData)
	{
		
		std::filesystem::create_directories(folder);

		std::string filename = folder + "/step_0_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream file(filename);

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

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	double lambda1 = (alpha + beta + sqrt((alpha - beta) * (alpha - beta) + 4 * gamma*gamma)) / 2.;
	double lambda2 = (alpha + beta - sqrt((alpha - beta) * (alpha - beta) + 4 * gamma*gamma)) / 2.;
	double w1 = 0.5 * (lambda1 / (lambda1 + lambda2));
	double w2 = 0.5 * (lambda2 / (lambda1 + lambda2));


	double b[8] = {
	alpha + 2 * gamma * (w1 - w2),   // E
	alpha + 2 * gamma * (w1 - w2),   // W
	beta + 2 * gamma * (w1 - w2),   // N
	beta + 2 * gamma * (w1 - w2),   // S
	2 * gamma * w2,                 // NE
	-2 * gamma * w1,                // NW
	-2 * gamma * w1,                // SE  
	2 * gamma * w2                  // SW  
	};

	const int d[8] = { +1, -1, +width, -width,
						+width + 1, +width - 1,
						-width + 1, -width - 1 };

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;
			if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			double bcenter = 0.0;
			for (int s = 0; s < 8; s++)
			{
				int q = k + d[s];
				M.insert(k, q) = -scale * b[s];
				bcenter += b[s];
			}
			M.insert(k, k) = 1.0 + scale * bcenter;
		}
	}

	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);

	int steps = static_cast<int>(std::round((Tend - Tstart) / tau));

	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		Eigen::VectorXd b(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2. + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 - h / 2. + j * h;

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

		// error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0 - h / 2 + j * h;
				int k = i * width + j;

				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;

				error_sum += diff * diff;
			}
		}
	}

	// EXPORT FINAL
	if (exportData)
	{
		std::string finalName = folder + "/final_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream finalFile(finalName);

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

double Image::S1_FBDS_Classic(bool exportData)
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

	//INIT 
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 - h / 2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 - h / 2 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}

	// prvotny error
	//for (int i = 0; i < height; i++)
	//{
	//	double y = -1.0 - h / 2 + i * h;
	//	for (int j = 0; j < width; j++)
	//	{
	//		double x = -1.0 - h / 2 + j * h;
	//		int k = i * width + j;

	//		double exact = u_exact(x, y, Tstart, &D[0][0]);
	//		double diff = u_old[k] - exact;
	//		error_sum += diff * diff;
	//	}
	//}

	std::string folder = "S1_FBDS_Classic";
	if (exportData)
	{
		std::filesystem::create_directories(folder);
		std::string filename = folder + "/step_0_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream file(filename);

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


	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	//smery E, W, N, S, NE, NW, SE, SW
	const int d[8] = { +1, -1, +width, -width,
			+width + 1, +width - 1,
			-width + 1, -width - 1 };

	//koeficienty pre jednotlive smery  v tom poradi ako su v poli d
	double b[8] = { alpha, alpha,
					beta, beta,
					gamma / 2.0, -gamma / 2.0,
					-gamma / 2.0, gamma / 2.0 };

	//tu ulozim tie koef forw/back nech to nemusim vsetko zvlast indexovat
	double b_forw[8], b_back[8];

	for (int i = 0; i < 8; i++)
	{
		b_forw[i] = std::max(0.0, b[i]);
		b_back[i] = std::min(0.0, b[i]);
	}

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;

			if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			double bcenter = 0.0;
			for (int s = 0; s < 8; s++)
			{
				int q = k + d[s];
				if (b_forw[s] > 0.0)
					M.insert(k, q) = -scale * b_forw[s];
				bcenter += b_forw[s];
			}
			M.insert(k, k) = 1.0 + scale * bcenter;
		}
	}
	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);
	if (solver.info() != Eigen::Success)
		std::cout << "decomposition failed" << std::endl;


	int steps = static_cast<int>(std::round((Tend - Tstart) / tau));

	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		std::vector<double> np_plus(size, 0.0);
		std::vector<double> np_minus(size, 0.0);

		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
					continue;

				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];
					double flux = b_back[s] * (u_old[k] - u_old[q]);
					if (u_old[k] > u_old[q])
						np_plus[k] += flux;
					else
						np_minus[k] += flux;
				}
				np_plus[k] = -np_plus[k];
				np_minus[k] = -np_minus[k];
			}
		}

		double umin = *std::min_element(u_old.begin(), u_old.end());
		double umax = *std::max_element(u_old.begin(), u_old.end());

		std::vector<double> theta_plus(size, 1.0);
		std::vector<double> theta_minus(size, 1.0);

		for (int k = 0; k < size; k++)
		{
			double up = u_old[k];
			//to checkuje ci to neni nula, lebo nechcem delit nulou
			theta_plus[k] = (np_plus[k] > 0.0) ? std::min(1.0, (umax - up) * h * h / (tau * np_plus[k])) : 1.0; 
			theta_minus[k] = (np_minus[k] < 0.0) ? std::min(1.0, (umin-up) * h * h / (tau * np_minus[k])) : 1.0;
		}



		Eigen::VectorXd rhs(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2. + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 - h / 2. + j * h;

				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
				{
					rhs(k) = u_exact(x, y, t_n + tau, &D[0][0]);
					continue;
				}

				rhs(k) = u_old[k];
				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];
					double theta_pq=1.0;
					if (u_old[k] > u_old[q])
						theta_pq = std::min(theta_plus[k], theta_minus[q]);
					else
						theta_pq = std::min(theta_minus[k], theta_plus[q]);

					double flux = theta_pq * b_back[s] * (u_old[k] - u_old[q]);
					rhs(k) -= scale * flux;
				}
			}
		}

		Eigen::VectorXd sol = solver.solve(rhs);
		if (solver.info() != Eigen::Success)
			std::cout << "solve failed at step " << n << std::endl;

		for (int k = 0; k < size; k++)
			u_old[k] = sol(k);



		// error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0 - h / 2 + j * h;
				int k = i * width + j;
				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;
				error_sum += diff * diff;
			}
		}
	}

	// EXPORT FINAL
	if (exportData)
	{
		std::string finalName = folder + "/final_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream finalFile(finalName);

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

double Image::S2_FBDS_Classic(bool exportData)
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

	//INIT 
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 - h / 2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 - h / 2 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}

	// prvotny error
	//for (int i = 0; i < height; i++)
	//{
	//	double y = -1.0 - h / 2 + i * h;
	//	for (int j = 0; j < width; j++)
	//	{
	//		double x = -1.0 - h / 2 + j * h;
	//		int k = i * width + j;

	//		double exact = u_exact(x, y, Tstart, &D[0][0]);
	//		double diff = u_old[k] - exact;
	//		error_sum += diff * diff;
	//	}
	//}

	std::string folder = "S2_FBDS_Classic";
	if (exportData)
	{
		std::filesystem::create_directories(folder);
		std::string filename = folder + "/step_0_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream file(filename);

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


	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	//smery E, W, N, S, NE, NW, SE, SW
	const int d[8] = { +1, -1, +width, -width,
			+width + 1, +width - 1,
			-width + 1, -width - 1 };

	//koeficienty pre jednotlive smery  v tom poradi ako su v poli d
	double b[8] = { alpha, alpha,
					beta, beta,
					gamma / 2.0, -gamma / 2.0,
					-gamma / 2.0, gamma / 2.0 };

	//tu ulozim tie koef forw/back nech to nemusim vsetko zvlast indexovat
	double b_forw[8], b_back[8];

	for (int i = 0; i < 8; i++)
	{
		b_forw[i] = std::max(0.0, b[i]);
		b_back[i] = std::min(0.0, b[i]);
	}

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;

			if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			double bcenter = 0.0;
			for (int s = 0; s < 8; s++)
			{
				int q = k + d[s];
				if (b_forw[s] > 0.0)
					M.insert(k, q) = -scale * b_forw[s];
				bcenter += b_forw[s];
			}
			M.insert(k, k) = 1.0 + scale * bcenter;
		}
	}
	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);
	if (solver.info() != Eigen::Success)
		std::cout << "decomposition failed" << std::endl;


	int steps = static_cast<int>(std::round((Tend - Tstart) / tau));

	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		double umin = std::numeric_limits<double>::max();
		double umax = std::numeric_limits<double>::lowest();

		for (int i = 1; i < height - 1; i++)
			for (int j = 1; j < width - 1; j++)
			{
				int k = i * width + j;
				umin = std::min(umin, u_old[k]);
				umax = std::max(umax, u_old[k]);
			}

		//prvy solve, theta_pq je vlastne =1, lebo nedavame limit
		Eigen::VectorXd rhs(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2. + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 - h / 2. + j * h;

				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
				{
					rhs(k) = u_exact(x, y, t_n + tau, &D[0][0]);
					continue;
				}

				rhs(k) = u_old[k];
				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];
					double flux = b_back[s] * (u_old[k] - u_old[q]); //(thetu ani nepisem lebo 1)
					rhs(k) -= scale * flux;
				}
			}
		}
		Eigen::VectorXd sol = solver.solve(rhs);//solve bez limitu

		//ulozim si ktore to porusuju
		std::vector<bool> violates(size, false);
		bool anyViolation = false;

		for (int i = 1; i < height - 1; i++)
		{
			for (int j = 1; j < width - 1; j++)
			{
				int k = i * width + j;
				if (sol(k) < umin || sol(k) > umax)
				{
					violates[k] = true;
					anyViolation = true;
				}
			}
		}

		//akoze teoreticky sa moze stat ze to sa nestane, right? idk
		int iter = 0;
		int maxIter = 20;

		while (anyViolation && iter < maxIter)
		{
			iter++;

			std::vector<bool> needsTheta(size, false);

			for (int k = 0; k < size; k++)
			{
				if (!violates[k]) continue;

				needsTheta[k] = true;

				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];

					if (q >= 0 && q < size)
					{
						int qi = q / width;
						int qj = q % width;

						if (qi > 0 && qi < height - 1 && qj > 0 && qj < width - 1)
							needsTheta[q] = true;
					}
				}
			}

			std::vector<double> np_plus(size, 0.0);
			std::vector<double> np_minus(size, 0.0);

			for (int i = 1; i < height - 1; i++)
			{
				for (int j = 1; j < width - 1; j++)
				{
					int k = i * width + j;

					if (!needsTheta[k]) continue;

					for (int s = 0; s < 8; s++)
					{
						int q = k + d[s];

						double flux = b_back[s] * (u_old[k] - u_old[q]);

						if (u_old[k] > u_old[q])
							np_plus[k] += flux;
						else
							np_minus[k] += flux;
					}

					np_plus[k] = -np_plus[k];
					np_minus[k] = -np_minus[k];
				}
			}

			std::vector<double> theta_plus(size, 1.0);
			std::vector<double> theta_minus(size, 1.0);

			for (int k = 0; k < size; k++)
			{
				if (!needsTheta[k]) continue;

				double up = u_old[k];

				theta_plus[k] = (np_plus[k] > 0.0)
					? std::min(1.0, (umax - up) * h * h / (tau * np_plus[k]))
					: 1.0;

				theta_minus[k] = (np_minus[k] < 0.0)
					? std::min(1.0, (umin - up) * h * h / (tau * np_minus[k]))
					: 1.0;
			}

			for (int i = 0; i < height; i++)
			{
				double y = -1.0 - h / 2.0 + i * h;

				for (int j = 0; j < width; j++)
				{
					int k = i * width + j;
					double x = -1.0 - h / 2.0 + j * h;

					if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
					{
						rhs(k) = u_exact(x, y, t_n + tau, &D[0][0]);
						continue;
					}

					rhs(k) = u_old[k];

					for (int s = 0; s < 8; s++)
					{
						int q = k + d[s];

						double theta_pq = 1.0;

						if (u_old[k] > u_old[q])
							theta_pq = std::min(theta_plus[k], theta_minus[q]);
						else
							theta_pq = std::min(theta_minus[k], theta_plus[q]);

						double flux = theta_pq * b_back[s] * (u_old[k] - u_old[q]);
						rhs(k) -= scale * flux;
					}
				}
			}

			sol = solver.solve(rhs);

			if (solver.info() != Eigen::Success)
				std::cout << "solve failed at step " << n << ", S2 iter " << iter << std::endl;

			// opatovny check na iteraciu
			std::fill(violates.begin(), violates.end(), false);
			anyViolation = false;

			for (int i = 1; i < height - 1; i++)
			{
				for (int j = 1; j < width - 1; j++)
				{
					int k = i * width + j;

					if (sol(k) < umin || sol(k) > umax)
					{
						violates[k] = true;
						anyViolation = true;
					}
				}
			}
		}

		//toto je rovnaku odtialto
		for (int k = 0; k < size; k++)
			u_old[k] = sol(k);

		// error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0 - h / 2 + j * h;
				int k = i * width + j;
				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;
				error_sum += diff * diff;
			}
		}
	}

	// EXPORT FINAL
	if (exportData)
	{
		std::string finalName = folder + "/final_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream finalFile(finalName);

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

double Image::S1_FBDS_ADCM(bool exportData)
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

	//INIT 
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 - h / 2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 - h / 2 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}


	std::string folder = "S1_FBDS_ADCM";
	if (exportData)
	{
		std::filesystem::create_directories(folder);
		std::string filename = folder + "/step_0_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream file(filename);

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


	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	//smery E, W, N, S, NE, NW, SE, SW
	const int d[8] = { +1, -1, +width, -width,
			+width + 1, +width - 1,
			-width + 1, -width - 1 };

	//koeficienty pre jednotlive smery  v tom poradi ako su v poli d
	double lambda1 = (alpha + beta + sqrt((alpha - beta) * (alpha - beta) + 4 * gamma * gamma)) / 2.;
	double lambda2 = (alpha + beta - sqrt((alpha - beta) * (alpha - beta) + 4 * gamma * gamma)) / 2.;
	double w1 = 0.5 * (lambda1 / (lambda1 + lambda2));
	double w2 = 0.5 * (lambda2 / (lambda1 + lambda2));


	double b[8] = {
	alpha + 2 * gamma * (w1 - w2),   // E
	alpha + 2 * gamma * (w1 - w2),   // W
	beta + 2 * gamma * (w1 - w2),   // N
	beta + 2 * gamma * (w1 - w2),   // S
	2 * gamma * w2,                 // NE
	-2 * gamma * w1,                // NW
	-2 * gamma * w1,                // SE  
	2 * gamma * w2                  // SW  
	};

	//tu ulozim tie koef forw/back nech to nemusim vsetko zvlast indexovat
	double b_forw[8], b_back[8];

	for (int i = 0; i < 8; i++)
	{
		b_forw[i] = std::max(0.0, b[i]);
		b_back[i] = std::min(0.0, b[i]);
	}

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;

			if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			double bcenter = 0.0;
			for (int s = 0; s < 8; s++)
			{
				int q = k + d[s];
				if (b_forw[s] > 0.0)
					M.insert(k, q) = -scale * b_forw[s];
				bcenter += b_forw[s];
			}
			M.insert(k, k) = 1.0 + scale * bcenter;
		}
	}
	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);
	if (solver.info() != Eigen::Success)
		std::cout << "decomposition failed" << std::endl;


	int steps = static_cast<int>(std::round((Tend - Tstart) / tau));

	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		std::vector<double> np_plus(size, 0.0);
		std::vector<double> np_minus(size, 0.0);

		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
					continue;

				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];
					double flux = b_back[s] * (u_old[k] - u_old[q]);
					if (u_old[k] > u_old[q])
						np_plus[k] += flux;
					else
						np_minus[k] += flux;
				}
				np_plus[k] = -np_plus[k];
				np_minus[k] = -np_minus[k];
			}
		}

		double umin = *std::min_element(u_old.begin(), u_old.end());
		double umax = *std::max_element(u_old.begin(), u_old.end());

		std::vector<double> theta_plus(size, 1.0);
		std::vector<double> theta_minus(size, 1.0);

		for (int k = 0; k < size; k++)
		{
			double up = u_old[k];
			//to checkuje ci to neni nula, lebo nechcem delit nulou
			theta_plus[k] = (np_plus[k] > 0.0) ? std::min(1.0, (umax - up) * h * h / (tau * np_plus[k])) : 1.0;
			theta_minus[k] = (np_minus[k] < 0.0) ? std::min(1.0, (umin - up) * h * h / (tau * np_minus[k])) : 1.0;
		}



		Eigen::VectorXd rhs(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2. + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 - h / 2. + j * h;

				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
				{
					rhs(k) = u_exact(x, y, t_n + tau, &D[0][0]);
					continue;
				}

				rhs(k) = u_old[k];
				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];
					double theta_pq = 1.0;
					if (u_old[k] > u_old[q])
						theta_pq = std::min(theta_plus[k], theta_minus[q]);
					else
						theta_pq = std::min(theta_minus[k], theta_plus[q]);

					double flux = theta_pq * b_back[s] * (u_old[k] - u_old[q]);
					rhs(k) -= scale * flux;
				}
			}
		}

		Eigen::VectorXd sol = solver.solve(rhs);
		if (solver.info() != Eigen::Success)
			std::cout << "solve failed at step " << n << std::endl;

		for (int k = 0; k < size; k++)
			u_old[k] = sol(k);



		// error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0 - h / 2 + j * h;
				int k = i * width + j;
				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;
				error_sum += diff * diff;
			}
		}
	}

	// EXPORT FINAL
	if (exportData)
	{
		std::string finalName = folder + "/final_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream finalFile(finalName);

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

double Image::S2_FBDS_ADCM(bool exportData)
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

	//INIT 
	for (int i = 0; i < height; i++)
	{
		double y = -1.0 - h / 2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 - h / 2 + j * h;
			u_old[i * width + j] = u_exact(x, y, Tstart, &D[0][0]);
		}
	}

	std::string folder = "S2_FBDS_ADCM";
	if (exportData)
	{
		std::filesystem::create_directories(folder);
		std::string filename = folder + "/step_0_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream file(filename);

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


	double scale = tau / (h * h);
	double alpha = D[0][0];
	double beta = D[1][1];
	double gamma = D[0][1];

	//smery E, W, N, S, NE, NW, SE, SW
	const int d[8] = { +1, -1, +width, -width,
			+width + 1, +width - 1,
			-width + 1, -width - 1 };

	//koeficienty pre jednotlive smery  v tom poradi ako su v poli d
	double lambda1 = (alpha + beta + sqrt((alpha - beta) * (alpha - beta) + 4 * gamma * gamma)) / 2.;
	double lambda2 = (alpha + beta - sqrt((alpha - beta) * (alpha - beta) + 4 * gamma * gamma)) / 2.;
	double w1 = 0.5 * (lambda1 / (lambda1 + lambda2));
	double w2 = 0.5 * (lambda2 / (lambda1 + lambda2));


	double b[8] = {
	alpha + 2 * gamma * (w1 - w2),   // E
	alpha + 2 * gamma * (w1 - w2),   // W
	beta + 2 * gamma * (w1 - w2),   // N
	beta + 2 * gamma * (w1 - w2),   // S
	2 * gamma * w2,                 // NE
	-2 * gamma * w1,                // NW
	-2 * gamma * w1,                // SE  
	2 * gamma * w2                  // SW  
	};

	//tu ulozim tie koef forw/back nech to nemusim vsetko zvlast indexovat
	double b_forw[8], b_back[8];

	for (int i = 0; i < 8; i++)
	{
		b_forw[i] = std::max(0.0, b[i]);
		b_back[i] = std::min(0.0, b[i]);
	}

	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int k = i * width + j;

			if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
			{
				M.insert(k, k) = 1.0;
				continue;
			}

			double bcenter = 0.0;
			for (int s = 0; s < 8; s++)
			{
				int q = k + d[s];
				if (b_forw[s] > 0.0)
					M.insert(k, q) = -scale * b_forw[s];
				bcenter += b_forw[s];
			}
			M.insert(k, k) = 1.0 + scale * bcenter;
		}
	}
	M.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);
	if (solver.info() != Eigen::Success)
		std::cout << "decomposition failed" << std::endl;


	int steps = static_cast<int>(std::round((Tend - Tstart) / tau));

	for (int n = 0; n < steps; n++)
	{
		double t_n = Tstart + n * tau;

		double umin = std::numeric_limits<double>::max();
		double umax = std::numeric_limits<double>::lowest();

		for (int i = 1; i < height - 1; i++)
		{
			for (int j = 1; j < width - 1; j++)
			{
				int k = i * width + j;
				umin = std::min(umin, u_old[k]);
				umax = std::max(umax, u_old[k]);
			}
		}

		//prvy solve, theta_pq je vlastne =1, lebo nedavame limit
		Eigen::VectorXd rhs(size);

		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2. + i * h;
			for (int j = 0; j < width; j++)
			{
				int k = i * width + j;
				double x = -1.0 - h / 2. + j * h;

				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
				{
					rhs(k) = u_exact(x, y, t_n + tau, &D[0][0]);
					continue;
				}

				rhs(k) = u_old[k];

				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];
					double flux = b_back[s] * (u_old[k] - u_old[q]);
					rhs(k) -= scale * flux;
				}
			}
		}

		Eigen::VectorXd sol = solver.solve(rhs);//solve bez limitu

		//ulozim si ktore to porusuju
		std::vector<bool> violates(size, false);
		bool anyViolation = false;

		double eps = 1e-12;

		for (int i = 1; i < height - 1; i++)
		{
			for (int j = 1; j < width - 1; j++)
			{
				int k = i * width + j;

				if (sol(k) < umin - eps || sol(k) > umax + eps)
				{
					violates[k] = true;
					anyViolation = true;
				}
			}
		}

		int iter = 0;
		int maxIter = 20;

		while (anyViolation && iter < maxIter)
		{
			iter++;

			std::vector<bool> needsTheta(size, false);

			for (int k = 0; k < size; k++)
			{
				if (!violates[k]) continue;

				needsTheta[k] = true;

				for (int s = 0; s < 8; s++)
				{
					int q = k + d[s];

					if (q >= 0 && q < size)
					{
						int qi = q / width;
						int qj = q % width;

						if (qi > 0 && qi < height - 1 && qj > 0 && qj < width - 1)
							needsTheta[q] = true;
					}
				}
			}

			std::vector<double> np_plus(size, 0.0);
			std::vector<double> np_minus(size, 0.0);

			for (int i = 1; i < height - 1; i++)
			{
				for (int j = 1; j < width - 1; j++)
				{
					int k = i * width + j;

					if (!needsTheta[k]) continue;

					for (int s = 0; s < 8; s++)
					{
						int q = k + d[s];

						double flux = b_back[s] * (u_old[k] - u_old[q]);

						if (u_old[k] > u_old[q])
							np_plus[k] += flux;
						else
							np_minus[k] += flux;
					}

					np_plus[k] = -np_plus[k];
					np_minus[k] = -np_minus[k];
				}
			}

			std::vector<double> theta_plus(size, 1.0);
			std::vector<double> theta_minus(size, 1.0);

			for (int k = 0; k < size; k++)
			{
				if (!needsTheta[k]) continue;

				double up = u_old[k];

				theta_plus[k] = (np_plus[k] > 0.0)
					? std::min(1.0, (umax - up) * h * h / (tau * np_plus[k]))
					: 1.0;

				theta_minus[k] = (np_minus[k] < 0.0)
					? std::min(1.0, (umin - up) * h * h / (tau * np_minus[k]))
					: 1.0;
			}

			for (int i = 0; i < height; i++)
			{
				double y = -1.0 - h / 2.0 + i * h;

				for (int j = 0; j < width; j++)
				{
					int k = i * width + j;
					double x = -1.0 - h / 2.0 + j * h;

					if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
					{
						rhs(k) = u_exact(x, y, t_n + tau, &D[0][0]);
						continue;
					}

					rhs(k) = u_old[k];

					for (int s = 0; s < 8; s++)
					{
						int q = k + d[s];

						double theta_pq = 1.0;

						if (u_old[k] > u_old[q])
							theta_pq = std::min(theta_plus[k], theta_minus[q]);
						else
							theta_pq = std::min(theta_minus[k], theta_plus[q]);

						double flux = theta_pq * b_back[s] * (u_old[k] - u_old[q]);
						rhs(k) -= scale * flux;
					}
				}
			}

			sol = solver.solve(rhs);

			if (solver.info() != Eigen::Success)
				std::cout << "solve failed at step " << n << ", S2 iter " << iter << std::endl;

			std::fill(violates.begin(), violates.end(), false);
			anyViolation = false;

			for (int i = 1; i < height - 1; i++)
			{
				for (int j = 1; j < width - 1; j++)
				{
					int k = i * width + j;

					if (sol(k) < umin - eps || sol(k) > umax + eps)
					{
						violates[k] = true;
						anyViolation = true;
					}
				}
			}
		}
		//toto je rovnaku odtialto
		for (int k = 0; k < size; k++)
			u_old[k] = sol(k);

		// error
		for (int i = 0; i < height; i++)
		{
			double y = -1.0 - h / 2 + i * h;
			for (int j = 0; j < width; j++)
			{
				double x = -1.0 - h / 2 + j * h;
				int k = i * width + j;
				double exact = u_exact(x, y, t_n + tau, &D[0][0]);
				double diff = u_old[k] - exact;
				error_sum += diff * diff;
			}
		}
	}

	// EXPORT FINAL
	if (exportData)
	{
		std::string finalName = folder + "/final_N" + std::to_string(N) +
			"_theta" + std::to_string(static_cast<int>(theta)) + ".csv";

		std::ofstream finalFile(finalName);

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
		D[0][0] = 0.489794;
		D[0][1] = -0.086213;
		D[1][0] = -0.086213;
		D[1][1] = 0.016206;
	}
	else if (theta == 22.5) {
		D[0][0] = 0.431191;
		D[0][1] = -0.178191;
		D[1][0] = -0.178191;
		D[1][1] = 0.074809;
	}
	else if (theta == 45.0) {
		D[0][0] = 0.253;
		D[0][1] = -0.252;
		D[1][0] = -0.252;
		D[1][1] = 0.253;
	}
	else if (theta == 5.0) {
		D[0][0] = 0.5;
		D[0][1] = -0.045;
		D[1][0] = -0.045;
		D[1][1] = 0.005;
	}
}

