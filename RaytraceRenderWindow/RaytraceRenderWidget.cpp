//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <random>
#include <QTimer>
// include the header file
#include "RaytraceRenderWidget.h"

#define N_THREADS 16
#define N_LOOPS 100
#define N_BOUNCES 5
#define TERMINATION_FACTOR 0.35f

// constructor
RaytraceRenderWidget::RaytraceRenderWidget
        (   
        // the geometric object to show
        std::vector<ThreeDModel>      *newTexturedObject,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // parent widget in visual hierarchy
        QWidget             *parent
        )
    // the : indicates variable instantiation rather than arbitrary code
    // it is considered good style to use it where possible
    : 
    // start by calling inherited constructor with parent widget's pointer
    QOpenGLWidget(parent),
    // then store the pointers that were passed in
    texturedObjects(newTexturedObject),
    renderParameters(newRenderParameters),
    raytraceScene(texturedObjects,renderParameters)
    { // constructor
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        restartRaytrace = false;
        QTimer *timer = new QTimer(this);
        connect(timer, &QTimer::timeout, this, &RaytraceRenderWidget::forceRepaint);
        timer->start(30);
        pthread_rwlock_init(&rw_lock, NULL);
    // leaves nothing to put into the constructor body
    } // constructor    
    // init lock

void RaytraceRenderWidget::forceRepaint(){
    update();
}
// destructor
RaytraceRenderWidget::~RaytraceRenderWidget()
    { // destructor
    // empty (for now)
    // all of our pointers are to data owned by another class
    // so we have no responsibility for destruction
    // and OpenGL cleanup is taken care of by Qt
    // destory lock
    pthread_rwlock_destroy(&rw_lock);
    } // destructor                                                                 

// called when OpenGL context is set up
void RaytraceRenderWidget::initializeGL()
    { // RaytraceRenderWidget::initializeGL()
	// this should remain empty
    } // RaytraceRenderWidget::initializeGL()

// called every time the widget is resized
void RaytraceRenderWidget::resizeGL(int w, int h)
    { // RaytraceRenderWidget::resizeGL()
    // resize the render image
    frameBuffer.Resize(w, h);
    } // RaytraceRenderWidget::resizeGL()
    
// called every time the widget needs painting
void RaytraceRenderWidget::paintGL()
    { // RaytraceRenderWidget::paintGL()
    // set background colour to white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    // and display the image
    glDrawPixels(frameBuffer.width, frameBuffer.height, GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer.block);
    } // RaytraceRenderWidget::paintGL()


    // routine that generates the image
void RaytraceRenderWidget::Raytrace()
{ // RaytraceRenderWidget::Raytrace()

    restartRaytrace = true;
    if(raytracingThread.joinable())
        raytracingThread.join();
    restartRaytrace = false;

    // thread safety:
    pthread_rwlock_wrlock(&rw_lock);
    //To make our lifes easier, lets calculate things on VCS.
    //So we need to process our scene to get a triangle soup in VCS.
    //IMPORTANT: You still need to complete the method that gets the modelview matrix in the scene class!
    raytraceScene.updateScene();

    //clear frame buffer before we start
    frameBuffer.clear(RGBAValue(0.0f, 0.0f, 0.0f,1.0f));
    floatbuffer = std::vector<std::vector<Homogeneous4>>(static_cast<unsigned int>(frameBuffer.height));
    for (std::vector<Homogeneous4> &row : floatbuffer) {
        row = std::vector<Homogeneous4>(static_cast<unsigned int>(frameBuffer.width));
    }
    pthread_rwlock_unlock(&rw_lock);

    raytracingThread= std::thread(&RaytraceRenderWidget::RaytraceThread,this);
    raytracingThread.detach();

} // RaytraceRenderWidget::Raytrace()
    

void RaytraceRenderWidget::RaytraceThread()
{
    int loops = renderParameters->monteCarloEnabled? N_LOOPS:1;
    std::cout << "I Will do " << loops << " loops" << std::endl;

    //Each pixel in parallel using openMP.
    for(int loop = 0; loop < loops; loop++){
        float aspect_radio = static_cast<float>(frameBuffer.width) /  static_cast<float>(frameBuffer.height);
        #pragma omp parallel for schedule(dynamic)
        for(int j = 0; j < frameBuffer.height; j++){
            for(int i = 0; i < frameBuffer.width; i++) {
                pthread_rwlock_rdlock(&rw_lock);
                //TODO: YOUR CODE GOES HERE
                Homogeneous4 color; // calculate your raytraced color here.

                // VCS, eye is at (0, 0, 0)
                // 1. get pixel's center coord into VCS, and build ray
                float dx = 0.5f, dy = 0.5f;
                if (renderParameters->monteCarloEnabled) {
                    dx = static_cast<float>(rand()) / RAND_MAX;
                    dy = static_cast<float>(rand()) / RAND_MAX;
                }
                float y = (((static_cast<float>(j) + dx) / static_cast<float>(frameBuffer.height)) - 0.5f) * 2.f;
                float x = (((static_cast<float>(i) + dy) / static_cast<float>(frameBuffer.width)) - 0.5f) * 2.f;
                // 2. set camera radio
                if (aspect_radio < 1.f) {
                    y *= aspect_radio;
                } else {
                    x *= aspect_radio;
                }
                Ray r(Cartesian3(0, 0, 0), Cartesian3(x, y, -1).unit(), Ray::Type::primary);

                // 3. get color from ray in the scene
                color = raytraceScene.colorFromRay(r, 1.f, 10);
                //Gamma correction....
                float gamma = 2.2f;
                //We already calculate everything in float, so we just do gamma correction before putting it integer format.
                float pre = float(loop)/float(loop+1), now = 1.f - pre;
                floatbuffer[j][i] = pre*floatbuffer[j][i] + now * color;
                color.x = pow(floatbuffer[j][i].x,1/gamma);
                color.y = pow(floatbuffer[j][i].y,1/gamma);
                color.z = pow(floatbuffer[j][i].z,1/gamma);
                frameBuffer[j][i] = RGBAValue(color.x*255.0f,color.y*255.0f,color.z*255.0f,255.0f);
                pthread_rwlock_unlock(&rw_lock);
            }
        }
        std::cout <<"[" << loop <<"]" << " Done a loop!" << std::endl;
        if(restartRaytrace){
            return;
        }
    }
    std::cout << "Done!" << std::endl;
}

// mouse-handling
void RaytraceRenderWidget::mousePressEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mousePressEvent()
    // store the button for future reference
    int whichButton = int(event->button());
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;

    
    // and we want to force mouse buttons to allow shift-click to be the same as right-click
    unsigned int modifiers = event->modifiers();
    
    // shift-click (any) counts as right click
    if (modifiers & Qt::ShiftModifier)
        whichButton = Qt::RightButton;
    
    // send signal to the controller for detailed processing
    emit BeginScaledDrag(whichButton, x,y);
    } // RaytraceRenderWidget::mousePressEvent()
    
void RaytraceRenderWidget::mouseMoveEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseMoveEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit ContinueScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseMoveEvent()
    
void RaytraceRenderWidget::mouseReleaseEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseReleaseEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit EndScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseReleaseEvent()
