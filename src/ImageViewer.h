#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets>
#include "ui_ImageViewer.h"
#include "ViewerWidget.h"
#include "Image.h"



class ImageViewer : public QMainWindow
{
	Q_OBJECT

public:
	ImageViewer(QWidget* parent = Q_NULLPTR);

private:
	Ui::ImageViewerClass* ui;
	ViewerWidget* vW;

	QButtonGroup* viewGroup;

	QSettings settings;
	QMessageBox msgBox;

	QImage original;
	Image img_proc;

	//ImageViewer Events
	void closeEvent(QCloseEvent* event);

	//Image functions
	bool openImage(QString filename);
	bool saveImage(QString filename);
	bool invertColors();
	bool showOriginal();
	bool showDamaged();
	bool showMask();

	void radioButtonSetup();

private slots:
	void on_actionOpen_triggered();
	void on_actionSave_as_triggered();
	void on_actionExit_triggered();
	void on_actionInvert_triggered();
	void onViewChanged(QAbstractButton* button);
	void on_pushButton_generateMask_clicked();
	void on_pushButton_laplace_clicked();
};
