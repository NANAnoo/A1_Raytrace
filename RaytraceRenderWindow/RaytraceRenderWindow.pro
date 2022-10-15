QT+=opengl
TEMPLATE = app
TARGET = RaytraceRenderWindow
INCLUDEPATH += .

# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

win32{
    LIBS += -lopengl32
}


#adding openMP
QMAKE_CXXFLAGS+= -fopenmp -Wall
LIBS += -fopenmp

# Input
HEADERS += ArcBall.h \
           ArcBallWidget.h \
           Cartesian3.h \
           Homogeneous4.h \
           Light.h \
           Material.h \
           Matrix4.h \
           Quaternion.h \
           Ray.h \
           RaytraceRenderWidget.h \
           RenderController.h \
           RenderParameters.h \
           RenderWidget.h \
           RenderWindow.h \
           RGBAImage.h \
           RGBAValue.h \
           Scene.h \
           ThreeDModel.h \
           Triangle.h
SOURCES += ArcBall.cpp \
           ArcBallWidget.cpp \
           Cartesian3.cpp \
           Homogeneous4.cpp \
           Light.cpp \
           Material.cpp \
           Ray.cpp \
           Renderparameters.cpp \
           Scene.cpp \
           ThreeDModel.cpp \
           Triangle.cpp \
           main.cpp \
           Matrix4.cpp \
           Quaternion.cpp \
           RaytraceRenderWidget.cpp \
           RenderController.cpp \
           RenderWidget.cpp \
           RenderWindow.cpp \
           RGBAImage.cpp \
           RGBAValue.cpp
