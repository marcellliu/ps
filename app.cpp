#include "app.h"
#include "UI.h"

app::app(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::app)
{
    ui->setupUi(this);

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

app::~app()
{
    delete ui;
}
