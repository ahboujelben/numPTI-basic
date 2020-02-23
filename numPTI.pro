#-------------------------------------------------
#
# Project created by QtCreator 2014-01-02T17:45:02
#
#-------------------------------------------------

QT       += core

QT       += gui

QT       += opengl

QT       += printsupport

TARGET = numPTI
CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++14

TEMPLATE = app

win32 {
    contains(QT_ARCH, i386) {
        #32 bit
        LIBS += -lopengl32 $$PWD/libs/glew/bin/Release/Win32/glew32.dll
        LIBS += -L$$PWD/libs/glew/lib/Release/Win32/ -lglew32
    } else {
        #64 bit
        LIBS += -lopengl32 $$PWD/libs/glew/bin/Release/x64/glew32.dll
        LIBS += -L$$PWD/libs/glew/lib/Release/x64/ -lglew32
    }
}

unix {
    LIBS += -lGLEW
}


SOURCES += main.cpp \
    pore.cpp \
    node.cpp \
    network.cpp \
    mainwindow.cpp \
    cluster.cpp \
    shader.cpp \
    worker.cpp \
    widget3d.cpp \
    hoshenKopelmann.cpp \
    loadData.cpp \
    element.cpp \
    solver.cpp \
    generationRegular.cpp \
    misc.cpp \
    block.cpp \
    particle.cpp \
    artificial.cpp \
    drugflow.cpp \
    particleflow.cpp \
    angiogenesis.cpp \
    angioFlow.cpp \
    parentvessel.cpp \
    retina.cpp \
    libs/qcustomplot/qcustomplot.cpp

HEADERS += \
    pore.h \
    node.h \
    network.h \
    tools.h \
    mainwindow.h \
    widget3d.h \
    cluster.h \
    worker.h \
    element.h \
    block.h \
    particle.h \
    libs/qcustomplot/qcustomplot.h

INCLUDEPATH += libs

FORMS += \
    mainwindow.ui

CONFIG += warn_off

RESOURCES += \
    resource.qrc
