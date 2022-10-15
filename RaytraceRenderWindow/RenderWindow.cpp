/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Window
//  -----------------------------
//
//  The render window class is really just a container
//  for tracking the visual hierarchy.  While it is
//  entirely possible to use Qt Creator, I try to avoid
//  over-commitment to it because I need to write code in
//  multiple environments, some of which are not well-suited
//  to IDEs in general, let alone Qt Creator
//
//  Also, by providing sample code, the didactic purpose of
//  showing how things fit together is better served.
//
/////////////////////////////////////////////////////////////////

#include "RenderWindow.h"
#include "RenderParameters.h"

// constructor
RenderWindow::RenderWindow
        (
        // the object to be rendered
        std::vector<ThreeDModel>          *newTexturedObject,
        // the model object storing render parameters
        RenderParameters        *newRenderParameters,
        // the title for the window (with default value)
        const char              *windowName
        )
    // call the inherited constructor
    // NULL indicates that this widget has no parent
    // i.e. that it is a top-level window
    :
    // member instantiation
    QWidget(NULL),
    texturedObject(newTexturedObject),
    renderParameters(newRenderParameters)
    { // RenderWindow::RenderWindow()
    // set the window's title
    setWindowTitle(QString(windowName));

    // initialise the grid layout
    windowLayout = new QGridLayout(this);

    // create all of the widgets, starting with the custom render widgets
    renderWidget                = new RenderWidget              (newTexturedObject,     newRenderParameters,        this);
    raytraceRenderWidget        = new RaytraceRenderWidget      (newTexturedObject,     newRenderParameters,        this);

    // construct custom arcball Widgets
    modelRotator                = new ArcBallWidget             (                       this);

    // construct standard QT widgets
    // check boxes
    fresnelRenderingBox     = new QCheckBox                 ("Fresnel",            this);
    interpolationBox     = new QCheckBox                 ("Interpolation",            this);
    phongshadingBox      = new QCheckBox                 ("Phong",            this);
    shadowBox            = new QCheckBox                 ("Shadow",            this);
    reflectionBox        = new QCheckBox                 ("Reflection",            this);
    refractionBox        = new QCheckBox                 ("Refraction",            this);
    monteCarloBox        = new QCheckBox                 ("Monte-Carlo",            this);


    // spatial sliders
    xTranslateSlider            = new QSlider                   (Qt::Horizontal,        this);
    secondXTranslateSlider      = new QSlider                   (Qt::Horizontal,        this);
    yTranslateSlider            = new QSlider                   (Qt::Vertical,          this);
    zTranslateSlider            = new QSlider                   (Qt::Vertical,          this);



    // labels for sliders and arcballs
    modelRotatorLabel           = new QLabel                    ("Model",               this);
    yTranslateLabel             = new QLabel                    ("Y",                   this);
    zoomLabel                   = new QLabel                    ("Zm",                  this);


    //raytracing things
    raytraceButton              = new QPushButton               ("Raytrace",            this);
    connect(raytraceButton,&QPushButton::released,
            this,&RenderWindow::handle_raytrace);
    // add all of the widgets to the grid               Row         Column      Row Span    Column Span

    // the top two widgets have to fit to the widgets stack between them
    int nStacked = 13;

    windowLayout->addWidget(renderWidget,               0,          1,          nStacked,   1           );
    windowLayout->addWidget(yTranslateSlider,           0,          2,          nStacked,   1           );
    windowLayout->addWidget(zTranslateSlider,                 0,          4,          nStacked,   1           );
    windowLayout->addWidget(raytraceRenderWidget,       0,          5,          nStacked,   1           );
    windowLayout->addWidget(raytraceButton,             0,          6,          nStacked,   1           );

    // the stack in the middle
    windowLayout->addWidget(modelRotator,               0,          3,          1,          1           );
    windowLayout->addWidget(modelRotatorLabel,          1,          3,          1,          1           );
    windowLayout->addWidget(interpolationBox,           2,         3,          1,          1           );
    windowLayout->addWidget(phongshadingBox,            3,         3,          1,          1           );
    windowLayout->addWidget(shadowBox,                  4,         3,          1,          1           );
    windowLayout->addWidget(reflectionBox,              5,         3,          1,          1           );
    windowLayout->addWidget(refractionBox,              6,         3,          1,          1           );
    windowLayout->addWidget(fresnelRenderingBox,       7,         3,          1,          1           );
    windowLayout->addWidget(monteCarloBox,              8,         3,          1,          1           );

    // Translate Slider Row
    windowLayout->addWidget(xTranslateSlider,           nStacked,   1,          1,          1           );
    windowLayout->addWidget(yTranslateLabel,            nStacked,   2,          1,          1           );
    // nothing in column 3
    windowLayout->addWidget(zoomLabel,                  nStacked,   4,          1,          1           );
    windowLayout->addWidget(secondXTranslateSlider,     nStacked,   5,          1,          1           );


    // now reset all of the control elements to match the render parameters passed in
    ResetInterface();
    } // RenderWindow::RenderWindow()

// routine to reset interface
// sets every visual control to match the model
// gets called by the controller after each change in the model
void RenderWindow::ResetInterface()
    { // RenderWindow::ResetInterface()
    // set check boxes
    fresnelRenderingBox    ->setChecked        (renderParameters   ->  fresnelRendering);
    interpolationBox    ->setChecked        (renderParameters   ->  interpolationRendering);
    shadowBox    ->setChecked        (renderParameters   ->  shadowsEnabled);
    phongshadingBox    ->setChecked        (renderParameters   ->  phongEnabled);
    reflectionBox    ->setChecked        (renderParameters   ->  reflectionEnabled);
    refractionBox    ->setChecked        (renderParameters   ->  refractionEnabled);
    monteCarloBox    ->setChecked        (renderParameters   ->  monteCarloEnabled);


    // set sliders
    // x & y translate are scaled to notional unit sphere in render widgets
    // but because the slider is defined as integer, we multiply by a 100 for all sliders
    xTranslateSlider        ->setMinimum        (int(TRANSLATE_MIN                               * PARAMETER_SCALING));
    xTranslateSlider        ->setMaximum        (int(TRANSLATE_MAX                               * PARAMETER_SCALING));
    xTranslateSlider        ->setValue          (int(renderParameters -> xTranslate              * PARAMETER_SCALING));

    secondXTranslateSlider  ->setMinimum        (int(TRANSLATE_MIN                               * PARAMETER_SCALING));
    secondXTranslateSlider  ->setMaximum        (int(TRANSLATE_MAX                               * PARAMETER_SCALING));
    secondXTranslateSlider  ->setValue          (int(renderParameters -> xTranslate              * PARAMETER_SCALING));

    yTranslateSlider        ->setMinimum        (int(TRANSLATE_MIN                               * PARAMETER_SCALING));
    yTranslateSlider        ->setMaximum        (int(TRANSLATE_MAX                               * PARAMETER_SCALING));
    yTranslateSlider        ->setValue          (int(renderParameters -> yTranslate              * PARAMETER_SCALING));

    // zoom slider is a logarithmic scale, so we want a narrow range
    zTranslateSlider  ->setMinimum        (int(TRANSLATE_MIN                               * PARAMETER_SCALING));
    zTranslateSlider  ->setMaximum        (int(TRANSLATE_MAX                               * PARAMETER_SCALING));
    zTranslateSlider  ->setValue          (int(renderParameters -> zTranslate              * PARAMETER_SCALING));

    // main lighting parameters are simple 0.0-1.0

    // now flag them all for update
    renderWidget            ->update();
    raytraceRenderWidget    ->update();
    modelRotator            ->update();
    xTranslateSlider        ->update();
    secondXTranslateSlider  ->update();
    yTranslateSlider        ->update();
    zTranslateSlider        ->update();
    fresnelRenderingBox     ->update();
    interpolationBox        ->update();
    phongshadingBox         ->update();
    interpolationBox        ->update();
    shadowBox               ->update();
    reflectionBox           ->update();
    refractionBox           ->update();
    monteCarloBox           ->update();

    } // RenderWindow::ResetInterface()

void RenderWindow::handle_raytrace()
{
     raytraceRenderWidget->Raytrace();
}
