TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    $$PWD/../freeimodpoly.cpp
HEADERS += \
    $$PWD/../freeimodpoly.h

INCLUDEPATH += $$PWD/../../../Vespucci/MinGW_libs/include
DEPENDPATH += $$PWD/../../../Vespucci/MinGW_libs/include

INCLUDEPATH += $$PWD/../FreeIModPoly
DEPENDPATH += $$PWD/../FreeIModPoly

#Armadillo
win32: LIBS += -L$$PWD/../../../Vespucci/MinGW_libs/lib/ -larmadillo
win32-g++: PRE_TARGETDEPS += $$PWD/../../../Vespucci/MinGW_libs/lib/libarmadillo.a

#HDF5
win32: LIBS += -L$$PWD/../../../Vespucci/MinGW_libs/lib/ -lhdf5
wind32-g++: PRE_TARGETDEPS += $$PWD/../../../Vespucci/MinGW_libs/lib/libhdf5.a

#ARPACK-NG
win32: LIBS += -L$$PWD/../../../Vespucci/MinGW_libs/lib/ -larpack
win32-g++: PRE_TARGETDEPS += $$PWD/../../../Vespucci/MinGW_libs/lib/libarpack.a

#OpenBLAS (linked dynamically because arpack links it dynamically)
win32: LIBS += -L$$PWD/../../../Vespucci/MinGW_libs/lib/ -llibopenblas
win32-g++: PRE_TARGETDEPS += $$PWD/../../../Vespucci/MinGW_libs/lib/libopenblas.dll.a

#Libgfortran
win32: LIBS += -L$$PWD/../MinGW_libs/lib/$$PWD/../../../Vespucci/MinGW_libs/lib/ -lgfortran
win32-g++: PRE_TARGETDEPS += $$PWD/../../../Vespucci/MinGW_libs/lib/libgfortran.a



