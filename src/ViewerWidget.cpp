#include   "ViewerWidget.h"

ViewerWidget::ViewerWidget(QSize imgSize, QWidget* parent)
	: QWidget(parent)
{
	setAttribute(Qt::WA_StaticContents);
	setMouseTracking(true);
	if (imgSize != QSize(0, 0)) {
		img = new QImage(imgSize, QImage::Format_ARGB32);
		img->fill(Qt::white);
		resizeWidget(img->size());
		setPainter();
		setDataPtr();
	}
}
ViewerWidget::~ViewerWidget()
{
	delete painter;
	delete img;
}
void ViewerWidget::resizeWidget(QSize size)
{
	this->resize(size);
	this->setMinimumSize(size);
	this->setMaximumSize(size);
}

//Image functions
bool ViewerWidget::setImage(const QImage& inputImg)
{
	if (img != nullptr) {
		delete painter;
		delete img;
	}
	img = new QImage(inputImg);
	if (!img) {
		return false;
	}
	resizeWidget(img->size());
	setPainter();
	setDataPtr();
	update();

	return true;
}
bool ViewerWidget::isEmpty()
{
	if (img == nullptr) {
		return true;
	}

	if (img->size() == QSize(0, 0)) {
		return true;
	}
	return false;
}

bool ViewerWidget::changeSize(int width, int height)
{
	QSize newSize(width, height);

	if (newSize != QSize(0, 0)) {
		if (img != nullptr) {
			delete painter;
			delete img;
		}

		img = new QImage(newSize, QImage::Format_ARGB32);
		if (!img) {
			return false;
		}
		img->fill(Qt::white);
		resizeWidget(img->size());
		setPainter();
		setDataPtr();
		update();
	}

	return true;
}

void ViewerWidget::setPixel(int x, int y, uchar r, uchar g, uchar b, uchar a)
{
	size_t startbyte = y * img->bytesPerLine() + x * 4;
	data[startbyte] = r;
	data[startbyte + 1] = g;
	data[startbyte + 2] = b;
	data[startbyte + 3] = a;
}
void ViewerWidget::setPixel(int x, int y, uchar val)
{
	if (val > 255) val = 255;
	if (val < 0) val = 0;
	data[y * img->bytesPerLine() + x] = val;
}
void ViewerWidget::setPixel(int x, int y, double val)
{
	if (val > 1) val = 1;
	if (val < 0) val = 0;
	setPixel(x, y, static_cast<uchar>(255 * val));
}
void ViewerWidget::setPixel(int x, int y, double valR, double valG, double valB, double valA)
{
	if (valR > 1) valR = 1;
	if (valG > 1) valG = 1;
	if (valB > 1) valB = 1;
	if (valA > 1) valA = 1;

	if (valR < 0) valR = 0;
	if (valG < 0) valG = 0;
	if (valB < 0) valB = 0;
	if (valA < 0) valA = 0;

	setPixel(x, y, static_cast<uchar>(255 * valR), static_cast<uchar>(255 * valG), static_cast<uchar>(255 * valB), static_cast<uchar>(255 * valA));
}

//Slots
void ViewerWidget::paintEvent(QPaintEvent* event)
{
	if (img != nullptr) {
		QPainter painter(this);
		QRect area = event->rect();
		painter.drawImage(area, *img, area);
	}
}