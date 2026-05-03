#include "Image.h"

void Image::computeD(double v1, double v2, double D[2][2])
{
	double norm2 = v1 * v1 + v2 * v2;

	if (norm2 < 1e-12)
	{
		D[0][0] = K1; D[0][1] = 0.0;
		D[1][0] = 0.0; D[1][1] = K2;
		return;
	}

	D[0][0] = (K1 * v1 * v1 + K2 * v2 * v2) / norm2; // alpha
	D[1][1] = (K1 * v2 * v2 + K2 * v1 * v1) / norm2; // beta
	D[0][1] = (K2 - K1) * v1 * v2 / norm2;       // gamma
	D[1][0] = D[0][1];
}

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

void Image::generateRandomImage(int n)
{
	width = n;
	height = n;
	size = width * height;

	imageData.resize(size);
	damaged.resize(size);
	laplace.resize(size);
	smoothed.resize(size);
	mask.resize(size);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	for (int i = 0; i < size; i++)
	{
		imageData[i] = dist(gen);
	}
}

void Image::vectorField(double x1, double x2, double& v1, double& v2, int fieldId)
{
	if (fieldId == 0)
	{
		v1 = 1.0;
		v2 = 1.0;
	}
	else if (fieldId == 1)
	{
		v1 = -x2;
		v2 = x1;
	}
	else
	{
		v1 = std::sin(M_PI * x2);
		v2 = std::cos(M_PI * x1);
	}
}

void Image::variableDCM()
{
	std::vector<double> original = imageData;

	imageData = original;
	variableDCM_forField(0, evolutionField1);

	imageData = original;
	variableDCM_forField(1, evolutionField2);

	imageData = original;
	variableDCM_forField(2, evolutionField3);

	imageData = evolutionField1.back(); 
}

void Image::variableDCM_forField(int fieldId, std::vector<std::vector<double>>& storage)
{
	if (width <= 0 || height <= 0 || imageData.empty())
		return;

	std::vector<double> old = imageData;
	std::vector<double> next = imageData;

	//sem pojdu vektorove polia
	std::vector<double> v1(size);
	std::vector<double> v2(size);
	for (int y = 0; y < height; y++)	//fill 'em
	{
		for (int x = 0; x < width; x++)
		{
			int id = y * width + x;

			double x1 = -1.0 + (x + 0.5) * (2.0 / width);
			double x2 = -1.0 + (y + 0.5) * (2.0 / height);

			vectorField(x1, x2, v1[id], v2[id], fieldId);
		}
	}


	storage.clear();
	storage.push_back(imageData);

	const double h = 1.0;
	const double scale = tau / (h * h);

	const int d[8] = {
		+1, -1, +width, -width,
		+width + 1, +width - 1,
		-width + 1, -width - 1
	};

	const int E = 0;
	const int W = 1;
	const int S = 2;
	const int Nn = 3;
	const int SE = 4;
	const int SW = 5;
	const int NE = 6;
	const int NW = 7;

	//build matrix
	Eigen::SparseMatrix<double, Eigen::RowMajor> M(size, size);
	M.reserve(Eigen::VectorXi::Constant(size, 9));

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int id = y * width + x;

			//boolovske pemenne na okraje, mohlo by to byt v ife normalne ale takto rovno vidim co kde je
			bool left = (x == 0);
			bool right = (x == width - 1);
			bool top = (y == 0);
			bool bottom = (y == height - 1);

			double vP1 = v1[id];
			double vP2 = v2[id];

			double alphaE = 0.0, gammaE = 0.0;
			double alphaW = 0.0, gammaW = 0.0;
			double betaN = 0.0, gammaN = 0.0;
			double betaS = 0.0, gammaS = 0.0;

			//toto sa len rataju tie vektorove polia na hrane , podla toho sa tvori ta matica
			if (!right)
			{
				double DE[2][2];
				//vektorove pole v strede bodu East
				double vE1 = v1[id + d[E]];
				double vE2 = v2[id + d[E]];
				//vyratam maticu D podla vektoroveho pola na HRANE teda priemer medzi P a E
				computeD(0.5 * (vP1 + vE1), 0.5 * (vP2 + vE2), DE);

				alphaE = DE[0][0];
				gammaE = DE[0][1];
			}

			if (!left)
			{
				double DW[2][2];
				double vW1 = v1[id + d[W]];
				double vW2 = v2[id + d[W]];

				computeD(0.5 * (vP1 + vW1), 0.5 * (vP2 + vW2), DW);

				alphaW = DW[0][0];
				gammaW = DW[0][1];
			}

			if (!top)
			{
				double DN[2][2];
				double vN1 = v1[id + d[Nn]];
				double vN2 = v2[id + d[Nn]];

				computeD(0.5 * (vP1 + vN1), 0.5 * (vP2 + vN2), DN);

				betaN = DN[1][1];
				gammaN = DN[0][1];
			}

			if (!bottom)
			{
				double DS[2][2];
				double vS1 = v1[id + d[S]];
				double vS2 = v2[id + d[S]];

				computeD(0.5 * (vP1 + vS1), 0.5 * (vP2 + vS2), DS);

				betaS = DS[1][1];
				gammaS = DS[0][1];
			}

			double bcenter = 0.0;

			// ---------------- VNUTRO ----------------
			if (!left && !right && !top && !bottom)
			{
				double bE = alphaE + (gammaN - gammaS) / 4.0;
				double bW = alphaW + (gammaS - gammaN) / 4.0;
				double bN = betaN + (gammaE - gammaW) / 4.0;
				double bS = betaS + (gammaW - gammaE) / 4.0;

				double bNE = (gammaE + gammaN) / 4.0;
				double bNW = -(gammaW + gammaN) / 4.0;
				double bSE = -(gammaE + gammaS) / 4.0;
				double bSW = (gammaW + gammaS) / 4.0;

				bcenter = bE + bW + bN + bS + bNE + bNW + bSE + bSW;

				M.insert(id, id + d[E]) = -scale * bE;
				M.insert(id, id + d[W]) = -scale * bW;
				M.insert(id, id + d[Nn]) = -scale * bN;
				M.insert(id, id + d[S]) = -scale * bS;
				M.insert(id, id + d[NE]) = -scale * bNE;
				M.insert(id, id + d[NW]) = -scale * bNW;
				M.insert(id, id + d[SE]) = -scale * bSE;
				M.insert(id, id + d[SW]) = -scale * bSW;
			}

			// ---------------- BOTTOM LEFT ----------------
			else if (bottom && left)
			{
				double bE = alphaE + (gammaN - gammaE) / 4.0;
				double bN = betaN + (gammaE - gammaN) / 4.0;
				double bNE = (gammaE + gammaN) / 4.0;

				bcenter = bE + bN + bNE;

				M.insert(id, id + d[E]) = -scale * bE;
				M.insert(id, id + d[Nn]) = -scale * bN;
				M.insert(id, id + d[NE]) = -scale * bNE;
			}

			// ---------------- BOTTOM RIGHT ----------------
			else if (bottom && right)
			{
				double bW = alphaW + (gammaW - gammaN) / 4.0;
				double bN = betaN + (gammaN - gammaW) / 4.0;
				double bNW = -(gammaW + gammaN) / 4.0;

				bcenter = bW + bN + bNW;

				M.insert(id, id + d[W]) = -scale * bW;
				M.insert(id, id + d[Nn]) = -scale * bN;
				M.insert(id, id + d[NW]) = -scale * bNW;
			}

			// ---------------- TOP LEFT ----------------
			else if (top && left)
			{
				double bE = alphaE + (gammaE - gammaS) / 4.0;
				double bS = betaS + (gammaS - gammaE) / 4.0;
				double bSE = -(gammaE + gammaS) / 4.0;

				bcenter = bE + bS + bSE;

				M.insert(id, id + d[E]) = -scale * bE;
				M.insert(id, id + d[S]) = -scale * bS;
				M.insert(id, id + d[SE]) = -scale * bSE;
			}

			// ---------------- TOP RIGHT ----------------
			else if (top && right)
			{
				double bW = alphaW + (gammaS - gammaW) / 4.0;
				double bS = betaS + (gammaW - gammaS) / 4.0;
				double bSW = (gammaW + gammaS) / 4.0;

				bcenter = bW + bS + bSW;

				M.insert(id, id + d[W]) = -scale * bW;
				M.insert(id, id + d[S]) = -scale * bS;
				M.insert(id, id + d[SW]) = -scale * bSW;
			}

			// ---------------- LEFT EDGE ----------------
			else if (left)
			{
				double bE = alphaE + (gammaN - gammaS) / 4.0;
				double bN = betaN + (gammaE - gammaN) / 4.0;
				double bS = betaS + (gammaS - gammaE) / 4.0;
				double bNE = (gammaE + gammaN) / 4.0;
				double bSE = -(gammaE + gammaS) / 4.0;

				bcenter = bE + bN + bS + bNE + bSE;

				M.insert(id, id + d[E]) = -scale * bE;
				M.insert(id, id + d[Nn]) = -scale * bN;
				M.insert(id, id + d[S]) = -scale * bS;
				M.insert(id, id + d[NE]) = -scale * bNE;
				M.insert(id, id + d[SE]) = -scale * bSE;
			}

			// ---------------- RIGHT EDGE ----------------
			else if (right)
			{
				double bW = alphaW + (gammaS - gammaN) / 4.0;
				double bN = betaN + (gammaN - gammaW) / 4.0;
				double bS = betaS + (gammaW - gammaS) / 4.0;
				double bNW = -(gammaW + gammaN) / 4.0;
				double bSW = (gammaW + gammaS) / 4.0;

				bcenter = bW + bN + bS + bNW + bSW;

				M.insert(id, id + d[W]) = -scale * bW;
				M.insert(id, id + d[Nn]) = -scale * bN;
				M.insert(id, id + d[S]) = -scale * bS;
				M.insert(id, id + d[NW]) = -scale * bNW;
				M.insert(id, id + d[SW]) = -scale * bSW;
			}

			// ---------------- TOP EDGE ----------------
			else if (top)
			{
				double bE = alphaE + (gammaE - gammaS) / 4.0;
				double bW = alphaW + (gammaS - gammaW) / 4.0;
				double bS = betaS + (gammaW - gammaE) / 4.0;
				double bSE = -(gammaE + gammaS) / 4.0;
				double bSW = (gammaW + gammaS) / 4.0;

				bcenter = bE + bW + bS + bSE + bSW;

				M.insert(id, id + d[E]) = -scale * bE;
				M.insert(id, id + d[W]) = -scale * bW;
				M.insert(id, id + d[S]) = -scale * bS;
				M.insert(id, id + d[SE]) = -scale * bSE;
				M.insert(id, id + d[SW]) = -scale * bSW;
			}

			// ---------------- BOTTOM EDGE ----------------
			else if (bottom)
			{
				double bE = alphaE + (gammaN - gammaE) / 4.0;
				double bW = alphaW + (gammaW - gammaN) / 4.0;
				double bN = betaN + (gammaE - gammaW) / 4.0;
				double bNE = (gammaE + gammaN) / 4.0;
				double bNW = -(gammaW + gammaN) / 4.0;

				bcenter = bE + bW + bN + bNE + bNW;

				M.insert(id, id + d[E]) = -scale * bE;
				M.insert(id, id + d[W]) = -scale * bW;
				M.insert(id, id + d[Nn]) = -scale * bN;
				M.insert(id, id + d[NE]) = -scale * bNE;
				M.insert(id, id + d[NW]) = -scale * bNW;
			}

			M.insert(id, id) = 1.0 + scale * bcenter;
		}
	}

	M.makeCompressed();

	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
	solver.compute(M);

	if (solver.info() != Eigen::Success) {
		std::cout << "Decomposition failed" << std::endl;
		return;
	}
	for (int step = 0; step < timeSteps; step++)
	{
		old = next;

		Eigen::VectorXd rhs(size);

		for (int i = 0; i < size; i++)
			rhs(i) = old[i];

		Eigen::VectorXd sol = solver.solve(rhs);

		if (solver.info() != Eigen::Success)
			return;

		for (int i = 0; i < size; i++)
		{
			next[i] = sol(i);

			if (next[i] < 0.0) next[i] = 0.0;
			if (next[i] > 1.0) next[i] = 1.0;
		}

		if ((step + 1) % 5 == 0)
			storage.push_back(next);
	}

	if (storage.empty() || storage.back() != next)
		storage.push_back(next);

}

double* Image::getEvolutionFrame(int fieldId, int step)
{
	if (fieldId == 0)
		return evolutionField1[step].data();

	if (fieldId == 1)
		return evolutionField2[step].data();

	return evolutionField3[step].data();
}

int Image::getEvolutionFrameCount(int fieldId) const
{
	if (fieldId == 0)
		return static_cast<int>(evolutionField1.size());

	if (fieldId == 1)
		return static_cast<int>(evolutionField2.size());

	return static_cast<int>(evolutionField3.size());
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
				else if (i == 0 && j >= 1 && j <= width - 2)	//dole
				{
					M.insert(id, id) = 4.;
					M.insert(id, id + width) = -2.0;
					M.insert(id, id - 1) = -1.0;
					M.insert(id, id + 1) = -1.0;
				}
				else if (i == height - 1 && j >= 1 && j <= width - 2) //hore
				{
					M.insert(id, id) = 4.;
					M.insert(id, id - width) = -2.0;
					M.insert(id, id - 1) = -1.0;
					M.insert(id, id + 1) = -1.0;
				}
				else if (j == 0 && i >= 1 && i <= height - 2)	//vlavo
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
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id + width) = -2.0;
				M.insert(id, id + 1) = -2.0;
			}
			else if (i == 0 && j == width - 1)	//vpravo dole
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id + width) = -2.0;
				M.insert(id, id - 1) = -2.0;
			}
			else if (i == height - 1 && j == 0)	//vlavo hore
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id - width) = -2.0;
				M.insert(id, id + 1) = -2.0;
			}
			else if (i == height - 1 && j == width - 1)	//vpravo hore
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id - width) = -2.0;
				M.insert(id, id - 1) = -2.0;
			}

			//HRANY
			else if (i == 0 && j >= 1 && j <= width - 2)	//dole
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id + width) = -2.0;
				M.insert(id, id - 1) = -1.0;
				M.insert(id, id + 1) = -1.0;
			}
			else if (i == height - 1 && j >= 1 && j <= width - 2) //hore
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id - width) = -2.0;
				M.insert(id, id - 1) = -1.0;
				M.insert(id, id + 1) = -1.0;
			}
			else if (j == 0 && i >= 1 && i <= height - 2)	//vlavo
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id - width) = -1.0;
				M.insert(id, id + width) = -1.0;
				M.insert(id, id + 1) = -2.0;
			}
			else if (j == width - 1 && i >= 1 && i <= height - 2)	//vpravo
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id - width) = -1.0;
				M.insert(id, id + width) = -1.0;
				M.insert(id, id - 1) = -2.0;
			}
			//ZVYSOK
			else
			{
				M.insert(id, id) = 4. + lambda;
				M.insert(id, id - width) = -1.0;
				M.insert(id, id + width) = -1.0;
				M.insert(id, id - 1) = -1.0;
				M.insert(id, id + 1) = -1.0;
			}
			b(id) = lambda * laplace[id];
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
		double y = -1.0 - h / 2 + i * h;
		for (int j = 0; j < width; j++)
		{
			double x = -1.0 - h / 2 + j * h;
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

