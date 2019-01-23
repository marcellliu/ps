#ifndef APP_H
#define APP_H

#include <QMainWindow>
#include "near_ps.h"

namespace Ui {
class app;
}

class app : public QMainWindow
{
    Q_OBJECT

public:
    explicit app(QWidget *parent = 0);
    ~app();

private:
    Ui::app *ui;
    near_ps *np;
};

#endif // APP_H
