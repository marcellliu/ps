#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //    QString filename = QFileDialog::getOpenFileName(this,QObject::tr("Open Mask Image"),".",QObject::tr("Image Files(*.png *.jpg *.raw)"));
    //    cv::Mat image = imread(filename.toLocal8Bit().data(),cv::IMREAD_GRAYSCALE);
    //    qDebug() << filename;
    QString dirname = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                           "/home",
                                                           QFileDialog::ShowDirsOnly
                                                           | QFileDialog::DontResolveSymlinks);
    qDebug() << dirname;


    QString filename = dirname+"/set.json";
    QFile jsonfile(filename.toLocal8Bit().data());
    if (jsonfile.exists())
        jsonfile.open(QIODevice::ReadOnly | QIODevice::Text);
    QByteArray data = jsonfile.readAll();
    jsonfile.close();

    QJsonParseError error;
    QJsonDocument jsonDocument = QJsonDocument::fromJson(data, &error);
    QJsonObject read_json = jsonDocument.object();


    np = new near_ps(read_json,10);
    np->init();
}

MainWindow::~MainWindow()
{
    delete ui;
}
