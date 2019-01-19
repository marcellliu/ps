#-------------------------------------------------
#
# Project created by QtCreator 2018-12-24T15:01:02
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = near_ps
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
        mainwindow.cpp \
    near_ps.cpp

HEADERS += \
        mainwindow.h \
    near_ps.h

FORMS += \
        mainwindow.ui

#INCLUDEPATH += /usr/include/opencv \
#               /user/include/opencv2

LIBS += -llapack -lblas -larmadillo -lopenblas
#LIBS += /usr/lib/x86_64-linux-gnu/libopencv_highgui.so \
#        /usr/lib/x86_64-linux-gnu/libopencv_core.so    \
#        /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so \
#        /usr/lib/x86_64-linux-gnu/libopencv_imgcodecs.so
