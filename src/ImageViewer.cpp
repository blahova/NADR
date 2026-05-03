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
	viewGroup->addButton(ui->radioButton_smooth);
	connect(viewGroup, QOverload<QAbstractButton*>::of(&QButtonGroup::buttonClicked), this,&ImageViewer::onViewChanged);

	fieldGroup = new QButtonGroup(this);
	fieldGroup->setExclusive(true);

	fieldGroup->addButton(ui->radioButton_field1, 0);
	fieldGroup->addButton(ui->radioButton_field2, 1);
	fieldGroup->addButton(ui->radioButton_field3, 2);
	ui->radioButton_field1->setChecked(true);
	ui->radioButton_field1->setEnabled(false);
	ui->radioButton_field2->setEnabled(false);
	ui->radioButton_field3->setEnabled(false);

	connect(fieldGroup, QOverload<int>::of(&QButtonGroup::idClicked), this, &ImageViewer::onFieldChanged);

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
		ui->pushButton_generateMask->setEnabled(true);

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

bool ImageViewer::showSmooth()
{
	if (vW->isEmpty())
		return false;

	double* procData = nullptr;

	procData = img_proc.getSmoothed();

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

bool ImageViewer::showData(double* procData)
{
	if (vW->isEmpty() || !procData)
		return false;

	int width = img_proc.getwidth();
	int height = img_proc.getheight();

	for (int y = 0; y < height; y++)
	{
		int row = y * width;

		for (int x = 0; x < width; x++)
		{
			int id = row + x;

			double val = procData[id];
			if (val < 0.0) val = 0.0;
			if (val > 1.0) val = 1.0;

			uchar pixelValue = static_cast<uchar>(val * 255.0 + 0.5);
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
	ui->pushButton_smooth->setEnabled(false);
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
	else if (button == ui->radioButton_smooth)
	{
		showSmooth();
	}
}

void ImageViewer::onFieldChanged(int id)
{
	currentField = id;

	int frames = img_proc.getEvolutionFrameCount(currentField);
	if (frames == 0)
		return;

	int frame = ui->verticalSlider_time->value();

	if (frame >= frames)
		frame = frames - 1;

	ui->verticalSlider_time->setMaximum(frames - 1);
	ui->verticalSlider_time->setValue(frame);

	showData(img_proc.getEvolutionFrame(currentField, frame));
}



// push buttony
void ImageViewer::on_pushButton_generateMask_clicked()
{
	ui->radioButton_mask->setEnabled(true);
	ui->radioButton_damaged->setEnabled(true);
	ui->radioButton_laplace->setEnabled(false);
	ui->radioButton_smooth->setEnabled(false);
	ui->pushButton_laplace->setEnabled(true);
	ui->pushButton_smooth->setEnabled(false);

	ui->radioButton_damaged->setChecked(true);

	img_proc.generateMask(ui->spinBox_maska->value());

	showDamaged();
}

void ImageViewer::on_pushButton_laplace_clicked()
{
	ui->radioButton_laplace->setEnabled(true);
	ui->radioButton_laplace->setChecked(true);
	ui->radioButton_smooth->setEnabled(false);
	ui->pushButton_smooth->setEnabled(true);

	img_proc.Laplace();

	showLaplace();
}

void ImageViewer::on_pushButton_smooth_clicked()
{
	ui->radioButton_smooth->setEnabled(true);
	ui->radioButton_smooth->setChecked(true);

	img_proc.Smooth(ui->doubleSpinBox_smooth->value());

	showSmooth();
}

void ImageViewer::on_pushButton_anisotropic_clicked()
{
	int metoda = ui->comboBox_method->currentIndex();

	int N = ui->comboBox_gridSize->currentText().toInt();
	double theta = ui->comboBox_theta->currentText().toDouble();

	img_proc.setN(N);
	img_proc.setTheta(theta);

	if (metoda == 0)
		img_proc.Anisotropic_Classic(true);
	else if(metoda==1)
		img_proc.Anisotropic_Modified(true);
	else if (metoda == 2)
		img_proc.S1_FBDS_Classic(true);
	else if (metoda == 3)
		img_proc.S2_FBDS_Classic(true);
	else if (metoda == 4)
		img_proc.S1_FBDS_ADCM(true);
	else if (metoda == 5)
		img_proc.S2_FBDS_ADCM(true);
	

}


void ImageViewer::on_pushButton_EOC_clicked()
{
	double theta = ui->comboBox_theta->currentText().toDouble();
	int metoda = ui->comboBox_method->currentIndex();

	std::vector<int> grids = { 20, 40, 80, 160 };
	std::vector<double> errors;

	std::cout << "ERRORS for theta " << theta << ":\n";

	for (int N : grids)
	{
		img_proc.setN(N);
		img_proc.setTheta(theta);

		double err;

		if (metoda == 0)
			err = img_proc.Anisotropic_Classic(false);
		else if (metoda == 1)
			err = img_proc.Anisotropic_Modified(false);
		else if (metoda == 2)
			err = img_proc.S1_FBDS_Classic(false);
		else if (metoda == 3)
			err = img_proc.S2_FBDS_Classic(false);
		else if (metoda == 4)
			err = img_proc.S1_FBDS_ADCM(false);
		else if (metoda == 5)
			err = img_proc.S2_FBDS_ADCM(false);
		errors.push_back(err);

		std::cout <<  err << std::endl;
	}

	std::cout << "\nEOC for theta "<<theta<<":\n";
	for (size_t i = 1; i < errors.size(); i++)
	{
		double eoc = log(errors[i - 1] / errors[i]) / log(2.0);
		std::cout <<  eoc << std::endl ;
	}
	std::cout <<std::endl;
}

//posledne zadanie
void ImageViewer::on_pushButton_run_clicked()
{
	double tau = ui->doubleSpinBox_tau->value();
	int steps = ui->spinBox_steps->value();
	int N = ui->spinBox_N->value();
	double K1 = ui->doubleSpinBox_K1->value();
	double K2 = ui->doubleSpinBox_K2->value();

	img_proc.generateRandomImage(N);
	vW->changeSize(N, N);
	showOriginal();

	img_proc.setTau(tau);
	img_proc.setTimeSteps(steps);
	img_proc.setK1(K1);
	img_proc.setK2(K2);

	img_proc.variableDCM();
	ui->radioButton_field1->setEnabled(true);
	ui->radioButton_field2->setEnabled(true);
	ui->radioButton_field3->setEnabled(true);

	int frames = img_proc.getEvolutionFrameCount(currentField);
	int lastFrame = frames - 1;

	ui->verticalSlider_time->setMinimum(0);
	ui->verticalSlider_time->setMaximum(lastFrame);
	ui->verticalSlider_time->setValue(lastFrame);
	ui->verticalSlider_time->setEnabled(true);

	showData(img_proc.getEvolutionFrame(currentField, lastFrame));

}

void ImageViewer::on_verticalSlider_time_valueChanged(int value)
{
	if (!ui->verticalSlider_time->isEnabled())
		return;

	if (img_proc.getEvolutionFrameCount(currentField) == 0)
		return;

	double* data = img_proc.getEvolutionFrame(currentField, value);

	showData(data);
}
