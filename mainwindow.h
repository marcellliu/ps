#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "near_ps.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    near_ps *np;
};

#endif // MAINWINDOW_H
