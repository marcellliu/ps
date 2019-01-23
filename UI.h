/********************************************************************************
** Form generated from reading UI file 'app.ui'
**
** Created by: Qt User Interface Compiler version 5.9.5
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_H
#define UI_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QOpenGLWidget>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_app
{
public:
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout_2;
    QWidget *menu;
    QVBoxLayout *verticalLayout;
    QToolButton *toolButton_4;
    QToolButton *toolButton;
    QToolButton *toolButton_3;
    QToolButton *toolButton_2;
    QSpacerItem *verticalSpacer;
    QSplitter *splitter;
    QScrollArea *scrollArea;
    QWidget *scrollAreaWidgetContents;
    QHBoxLayout *horizontalLayout_3;
    QStackedWidget *Display_win;
    QWidget *page;
    QVBoxLayout *verticalLayout_3;
    QOpenGLWidget *openGLWidget;
    QWidget *page_2;
    QWidget *Cout_Win;
    QVBoxLayout *verticalLayout_2;
    QTextBrowser *textBrowser;
    QWidget *Bottom_Menu;
    QHBoxLayout *horizontalLayout;
    QToolButton *toolButton_7;
    QToolButton *toolButton_5;
    QToolButton *toolButton_6;
    QSpacerItem *horizontalSpacer;
    QToolButton *Text_Win;
    QMenuBar *menubar;

    void setupUi(QMainWindow *app)
    {
        if (app->objectName().isEmpty())
            app->setObjectName(QStringLiteral("app"));
        app->resize(800, 600);
        centralwidget = new QWidget(app);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        horizontalLayout_2 = new QHBoxLayout(centralwidget);
        horizontalLayout_2->setSpacing(0);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        menu = new QWidget(centralwidget);
        menu->setObjectName(QStringLiteral("menu"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(menu->sizePolicy().hasHeightForWidth());
        menu->setSizePolicy(sizePolicy);
        verticalLayout = new QVBoxLayout(menu);
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        toolButton_4 = new QToolButton(menu);
        toolButton_4->setObjectName(QStringLiteral("toolButton_4"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(toolButton_4->sizePolicy().hasHeightForWidth());
        toolButton_4->setSizePolicy(sizePolicy1);

        verticalLayout->addWidget(toolButton_4);

        toolButton = new QToolButton(menu);
        toolButton->setObjectName(QStringLiteral("toolButton"));
        sizePolicy1.setHeightForWidth(toolButton->sizePolicy().hasHeightForWidth());
        toolButton->setSizePolicy(sizePolicy1);

        verticalLayout->addWidget(toolButton);

        toolButton_3 = new QToolButton(menu);
        toolButton_3->setObjectName(QStringLiteral("toolButton_3"));
        sizePolicy1.setHeightForWidth(toolButton_3->sizePolicy().hasHeightForWidth());
        toolButton_3->setSizePolicy(sizePolicy1);

        verticalLayout->addWidget(toolButton_3);

        toolButton_2 = new QToolButton(menu);
        toolButton_2->setObjectName(QStringLiteral("toolButton_2"));
        sizePolicy1.setHeightForWidth(toolButton_2->sizePolicy().hasHeightForWidth());
        toolButton_2->setSizePolicy(sizePolicy1);

        verticalLayout->addWidget(toolButton_2);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        horizontalLayout_2->addWidget(menu);

        splitter = new QSplitter(centralwidget);
        splitter->setObjectName(QStringLiteral("splitter"));
        splitter->setMidLineWidth(0);
        splitter->setOrientation(Qt::Vertical);
        splitter->setChildrenCollapsible(false);
        scrollArea = new QScrollArea(splitter);
        scrollArea->setObjectName(QStringLiteral("scrollArea"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(99);
        sizePolicy2.setHeightForWidth(scrollArea->sizePolicy().hasHeightForWidth());
        scrollArea->setSizePolicy(sizePolicy2);
        scrollArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QStringLiteral("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 760, 600));
        QSizePolicy sizePolicy3(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(scrollAreaWidgetContents->sizePolicy().hasHeightForWidth());
        scrollAreaWidgetContents->setSizePolicy(sizePolicy3);
        horizontalLayout_3 = new QHBoxLayout(scrollAreaWidgetContents);
        horizontalLayout_3->setSpacing(0);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        Display_win = new QStackedWidget(scrollAreaWidgetContents);
        Display_win->setObjectName(QStringLiteral("Display_win"));
        sizePolicy3.setHeightForWidth(Display_win->sizePolicy().hasHeightForWidth());
        Display_win->setSizePolicy(sizePolicy3);
        page = new QWidget();
        page->setObjectName(QStringLiteral("page"));
        sizePolicy3.setHeightForWidth(page->sizePolicy().hasHeightForWidth());
        page->setSizePolicy(sizePolicy3);
        verticalLayout_3 = new QVBoxLayout(page);
        verticalLayout_3->setSpacing(0);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        openGLWidget = new QOpenGLWidget(page);
        openGLWidget->setObjectName(QStringLiteral("openGLWidget"));
        sizePolicy3.setHeightForWidth(openGLWidget->sizePolicy().hasHeightForWidth());
        openGLWidget->setSizePolicy(sizePolicy3);
        openGLWidget->setMinimumSize(QSize(600, 600));

        verticalLayout_3->addWidget(openGLWidget);

        Display_win->addWidget(page);
        page_2 = new QWidget();
        page_2->setObjectName(QStringLiteral("page_2"));
        Display_win->addWidget(page_2);

        horizontalLayout_3->addWidget(Display_win);

        scrollArea->setWidget(scrollAreaWidgetContents);
        splitter->addWidget(scrollArea);
        Cout_Win = new QWidget(splitter);
        Cout_Win->setObjectName(QStringLiteral("Cout_Win"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(1);
        sizePolicy4.setHeightForWidth(Cout_Win->sizePolicy().hasHeightForWidth());
        Cout_Win->setSizePolicy(sizePolicy4);
        verticalLayout_2 = new QVBoxLayout(Cout_Win);
        verticalLayout_2->setSpacing(0);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        textBrowser = new QTextBrowser(Cout_Win);
        textBrowser->setObjectName(QStringLiteral("textBrowser"));
        QSizePolicy sizePolicy5(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(textBrowser->sizePolicy().hasHeightForWidth());
        textBrowser->setSizePolicy(sizePolicy5);

        verticalLayout_2->addWidget(textBrowser);

        Bottom_Menu = new QWidget(Cout_Win);
        Bottom_Menu->setObjectName(QStringLiteral("Bottom_Menu"));
        QSizePolicy sizePolicy6(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(Bottom_Menu->sizePolicy().hasHeightForWidth());
        Bottom_Menu->setSizePolicy(sizePolicy6);
        horizontalLayout = new QHBoxLayout(Bottom_Menu);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        toolButton_7 = new QToolButton(Bottom_Menu);
        toolButton_7->setObjectName(QStringLiteral("toolButton_7"));

        horizontalLayout->addWidget(toolButton_7);

        toolButton_5 = new QToolButton(Bottom_Menu);
        toolButton_5->setObjectName(QStringLiteral("toolButton_5"));

        horizontalLayout->addWidget(toolButton_5);

        toolButton_6 = new QToolButton(Bottom_Menu);
        toolButton_6->setObjectName(QStringLiteral("toolButton_6"));

        horizontalLayout->addWidget(toolButton_6);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        Text_Win = new QToolButton(Bottom_Menu);
        Text_Win->setObjectName(QStringLiteral("Text_Win"));

        horizontalLayout->addWidget(Text_Win);


        verticalLayout_2->addWidget(Bottom_Menu);

        splitter->addWidget(Cout_Win);

        horizontalLayout_2->addWidget(splitter);

        app->setCentralWidget(centralwidget);
        menubar = new QMenuBar(app);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 800, 28));
        app->setMenuBar(menubar);

        retranslateUi(app);

        QMetaObject::connectSlotsByName(app);
    } // setupUi

    void retranslateUi(QMainWindow *app)
    {
        app->setWindowTitle(QApplication::translate("app", "MainWindow", Q_NULLPTR));
        toolButton_4->setText(QString());
        toolButton->setText(QString());
        toolButton_3->setText(QString());
        toolButton_2->setText(QString());
        toolButton_7->setText(QString());
        toolButton_5->setText(QApplication::translate("app", "Status Info.", Q_NULLPTR));
        toolButton_6->setText(QString());
        Text_Win->setText(QApplication::translate("app", "...", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class app: public Ui_app {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_H
