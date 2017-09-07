#-------------------------------------------------
#
# Project created by QtCreator 2017-08-25T16:19:36
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = funcConnGUI
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    editlist.cpp \
    editfeat.cpp \
    choosecombinations.cpp

HEADERS  += mainwindow.h \
    editlist.h \
    editfeat.h \
    choosecombinations.h

FORMS    += mainwindow.ui \
    editlist.ui \
    editfeat.ui \
    choosecombinations.ui

DISTFILES +=
