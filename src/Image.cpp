#include "Image.h"


Image::Image(uchar* data, int w, int h, int bytesPerLine)
{
	width = w;
	height = h;
	size = w * h;

	imageData.resize(size);
	damaged.resize(size);
	mask.resize(size);
	laplace.resize(size);

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
	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

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
