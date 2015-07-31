TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    $$PWD/../freeimodpoly.cpp
HEADERS += \
    $$PWD/../freeimodpoly.h

INCLUDEPATH += $$PWD/../
DEPENDPATH += $$PWD/../

INCLUDEPATH += $$PWD/../../../../../usr/include
DEPENDPATH += $$PWD/../../../../../usr/include



#Armadillo
win32: LIBS += -L$$PWD/../../../../../usr/lib/ -larmadillo
win32-g++: PRE_TARGETDEPS += $$PWD/../../../../../usr/lib/libarmadillo.a

#HDF5
win32: LIBS += -L$$PWD/../../../../../usr/lib/ -lhdf5
wind32-g++: PRE_TARGETDEPS += $$PWD/../../../../../usr/lib/libhdf5.a

#ARPACK-NG
win32: LIBS += -L$$PWD/../../../../../usr/lib/ -larpack
win32-g++: PRE_TARGETDEPS += $$PWD/../../../../../usr/lib/libarpack.a

#OpenBLAS (linked dynamically because arpack links it dynamically)
win32: LIBS += -L$$PWD/../../../../../usr/lib/ -llibopenblas
win32-g++: PRE_TARGETDEPS += $$PWD/../../../../../usr/lib/libopenblas.dll.a

#Libgfortran
win32: LIBS += -L$$PWD/../../../../../usr/lib/ -lgfortran
win32-g++: PRE_TARGETDEPS += $$PWD/../../../../../usr/lib/libgfortran.a

