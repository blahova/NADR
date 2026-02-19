#include "ImageViewer.h"



ImageViewer::ImageViewer(QWidget* parent)
	: QMainWindow(parent), ui(new Ui::ImageViewerClass)
{
	ui->setupUi(this);
	vW = new ViewerWidget(QSize(500, 500));
	ui->scrollArea->setWidget(vW);

	ui->scrollArea->setBackgroundRole(QPalette::Dark);
	ui->scrollArea->setWidgetResizable(true);
	ui->scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	ui->scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);

	vW->setObjectName("ViewerWidget");

	viewGroup = new QButtonGroup(this);
	viewGroup->setExclusive(true);

	viewGroup->addButton(ui->radioButton_original);
	viewGroup->addButton(ui->radioButton_damaged);
	viewGroup->addButton(ui->radioButton_mask);
	viewGroup->addButton(ui->radioButton_laplace);
	connect(viewGroup, QOverload<QAbstractButton*>::of(&QButtonGroup::buttonClicked), this,&ImageViewer::onViewChanged);

	radioButtonSetup();
}

//ImageViewer Events
void ImageViewer::closeEvent(QCloseEvent* event)
{
	if (QMessageBox::Yes == QMessageBox::question(this, "Close Confirmation", "Are you sure you want to exit?", QMessageBox::Yes | QMessageBox::No))
	{
		event->accept();
	}
	else {
		event->ignore();
	}
}

//Image functions
bool ImageViewer::openImage(QString filename)
{
	QImage loadedImg(filename);
	if (!loadedImg.isNull()) {
		img_proc = Image(loadedImg.bits(), loadedImg.width(), loadedImg.height(), loadedImg.bytesPerLine());

		radioButtonSetup();

		qDebug() << loadedImg.format()
			<< loadedImg.width()
			<< loadedImg.height()
			<< loadedImg.bytesPerLine();


		return vW->setImage(loadedImg);
	}
	return false;
}
bool ImageViewer::saveImage(QString filename)
{
	QFileInfo fi(filename);
	QString extension = fi.completeSuffix();

	QImage* img = vW->getImage();
	return img->save(filename, extension.toStdString().c_str());
}

bool ImageViewer::invertColors()
{
	if (vW->isEmpty()) {
		return false;
	}

	uchar* data = vW->getData();

	int row = vW->getImage()->bytesPerLine();
	int depth = vW->getImage()->depth();

	for (int i = 0; i < vW->getImgHeight(); i++)
	{
		for (int j = 0; j < vW->getImgWidth(); j++)
		{
			//Grayscale
			if (depth == 8) {
				vW->setPixel(j, i, static_cast<uchar>(255 - data[i * row + j]));
			}
			//RGBA
			else {
				uchar r = static_cast<uchar>(255 - data[i * row + j * 4]);
				uchar g = static_cast<uchar>(255 - data[i * row + j * 4 + 1]);
				uchar b = static_cast<uchar>(255 - data[i * row + j * 4 + 2]);
				vW->setPixel(j, i, r, g, b);
			}
		}
	}
	vW->update();
	return true;
}


//show functions
bool ImageViewer::showOriginal()
{
	if (vW->isEmpty())
		return false;

	double* procData = nullptr;

	procData = img_proc.getImageData();

	if (!procData)
		return false;

	int width = img_proc.getwidth();
	int height = img_proc.getheight();

	for (int y = 0; y < height; y++) {
		int row = y * width;
		for (int x = 0; x < width; x++) {
			int id = row + x;

			uchar pixelValue = static_cast<uchar>(procData[id] * 255.0 + 0.5);

			vW->setPixel(x, y, pixelValue);
		}
	}

	vW->update();
	return true;
}

bool ImageViewer::showDamaged()
{
	if (vW->isEmpty())
		return false;

	double* procData = nullptr;

	procData = img_proc.getDamaged();

	if (!procData)
		return false;

	int width = img_proc.getwidth();
	int height = img_proc.getheight();

	for (int y = 0; y < height; y++) {
		int row = y * width;
		for (int x = 0; x < width; x++) {
			int id = row + x;

			uchar pixelValue =static_cast<uchar>(procData[id] * 255.0 + 0.5);

			vW->setPixel(x, y, pixelValue);
		}
	}

	vW->update();
	return true;
}

bool ImageViewer::showMask()
{
	if (vW->isEmpty())
		return false;

	uint8_t* procData = nullptr;

	procData = img_proc.getMask();

	if (!procData)
		return false;

	int width = img_proc.getwidth();
	int height = img_proc.getheight();

	for (int y = 0; y < height; y++) {
		int row = y * width;
		for (int x = 0; x < width; x++) {
			int id = row + x;

			uchar pixelValue = static_cast<uchar>(procData[id] * 255.0 + 0.5);

			vW->setPixel(x, y, pixelValue);
		}
	}

	vW->update();
	return true;
}

bool ImageViewer::showLaplace()
{
	if (vW->isEmpty())
		return false;

	double* procData = nullptr;

	procData = img_proc.getLaplace();

	if (!procData)
		return false;

	int width = img_proc.getwidth();
	int height = img_proc.getheight();

	for (int y = 0; y < height; y++) {
		int row = y * width;
		for (int x = 0; x < width; x++) {
			int id = row + x;

			uchar pixelValue = static_cast<uchar>(procData[id] * 255.0 + 0.5);

			vW->setPixel(x, y, pixelValue);
		}
	}

	vW->update();
	return true;
}


void ImageViewer::radioButtonSetup()
{
	ui->radioButton_mask->setEnabled(false);
	ui->radioButton_damaged->setEnabled(false);
	ui->radioButton_laplace->setEnabled(false);
	ui->radioButton_original->setChecked(true);

	ui->pushButton_laplace->setEnabled(false);
}


//Slots
void ImageViewer::on_actionOpen_triggered()
{
	QString folder = settings.value("folder_img_load_path", "").toString();

	QString fileFilter = "Image data (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm .*xbm .* xpm);;All files (*)";
	QString fileName = QFileDialog::getOpenFileName(this, "Load image", folder, fileFilter);
	if (fileName.isEmpty()) { return; }

	QFileInfo fi(fileName);
	settings.setValue("folder_img_load_path", fi.absoluteDir().absolutePath());

	if (!openImage(fileName)) {
		msgBox.setText("Unable to open image.");
		msgBox.setIcon(QMessageBox::Warning);
		msgBox.exec();
	}
}
void ImageViewer::on_actionSave_as_triggered()
{
	QString folder = settings.value("folder_img_save_path", "").toString();

	QString fileFilter = "Image data (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm .*xbm .* xpm);;All files (*)";
	QString fileName = QFileDialog::getSaveFileName(this, "Save image", folder, fileFilter);
	if (!fileName.isEmpty()) {
		QFileInfo fi(fileName);
		settings.setValue("folder_img_save_path", fi.absoluteDir().absolutePath());

		if (!saveImage(fileName)) {
			msgBox.setText("Unable to save image.");
			msgBox.setIcon(QMessageBox::Warning);
		}
		else {
			msgBox.setText(QString("File %1 saved.").arg(fileName));
			msgBox.setIcon(QMessageBox::Information);
		}
		msgBox.exec();
	}
}
void ImageViewer::on_actionExit_triggered()
{
	this->close();
}

void ImageViewer::on_actionInvert_triggered()
{
	invertColors();
}

void ImageViewer::onViewChanged(QAbstractButton* button)
{
	if (button == ui->radioButton_original)
	{
		showOriginal();
	}
	else if (button == ui->radioButton_damaged)
	{
		showDamaged();
	}
	else if (button == ui->radioButton_mask)
	{
		showMask();
	}
	else if (button == ui->radioButton_laplace)
	{
		showLaplace();
	}
}



// push buttony
void ImageViewer::on_pushButton_generateMask_clicked()
{
	ui->radioButton_mask->setEnabled(true);
	ui->radioButton_damaged->setEnabled(true);
	ui->radioButton_laplace->setEnabled(false);
	ui->pushButton_laplace->setEnabled(true);

	ui->radioButton_damaged->setChecked(true);

	img_proc.generateMask(ui->spinBox_maska->value());

	showDamaged();
}

void ImageViewer::on_pushButton_laplace_clicked()
{
	ui->radioButton_laplace->setEnabled(true);
	ui->radioButton_laplace->setChecked(true);

	img_proc.Laplace();

	showLaplace();
}
