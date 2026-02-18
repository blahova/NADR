#include "Image.h"


Image::Image(uchar* data, int w, int h)
{
	width = w;
	height = h;
	size = h * w;

	imageData.resize(size);
	damaged.resize(size);
	mask.resize(size);
	laplace.resize(size);

	for (int i = 0; i < height; i++) {
		int row = i * width;
		for (int j = 0; j < width; j++) {
			int id = row + j;
			imageData[id] = static_cast<double>(data[id]) / 255.0;
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
