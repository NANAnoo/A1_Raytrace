/////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Render Parameters
//  -----------------------------
//
//  This is part of the "model" in the MVC paradigm
//  We separate out the render parameters from the object rendered
//
/////////////////////////////////////////////////////////////////

// include guard
#ifndef RENDER_PARAMETERS_H
#define RENDER_PARAMETERS_H

#include "Matrix4.h"
#include "Light.h"
#include <vector>

//here not to break the includes
class ThreeDModel;
#include "ThreeDModel.h"

// class for the render parameters
class RenderParameters
    { // class RenderParameters
    public:

    // we have a widget with an arcball which stores the rotation
    // we'll be lazy and leave it there instead of moving it here

    // we store x & y translations
    float xTranslate, yTranslate,zTranslate;



    std::vector<Light*> lights;

    // we will want two homogeneous matrices holding the rotations for the
    // model and the light.
    Matrix4 rotationMatrix;

    // and the booleans
    bool interpolationRendering;
    bool phongEnabled;
    bool fresnelRendering;
    bool shadowsEnabled;
    bool reflectionEnabled;
    bool refractionEnabled;
    bool monteCarloEnabled;

    bool centreObject;

    bool orthoProjection;


    // constructor
    RenderParameters()
        :
        xTranslate(0.0),
        yTranslate(0.0),
        zTranslate(0.0),
        interpolationRendering(false),
        phongEnabled(false),
        fresnelRendering(false),
        shadowsEnabled(false),
        reflectionEnabled(false),
        refractionEnabled(false),
        monteCarloEnabled(false),
        centreObject(false),
        orthoProjection(false)
        { // constructor

        // because we are paranoid, we will initialise the matrices to the identity
        rotationMatrix.SetIdentity();
        } // constructor

    // accessor for scaledXTranslate

    ~RenderParameters(){
        for (unsigned int i = 0;i<lights.size();i++) {
            delete(lights.at(i));
        }
    }

    void findLights(std::vector<ThreeDModel> objects);

    }; // class RenderParameters

// now define some macros for bounds on parameters
#define TRANSLATE_MIN -1.0f
#define TRANSLATE_MAX 1.0f

#define ZOOM_SCALE_LOG_MIN -2.0f
#define ZOOM_SCALE_LOG_MAX 2.0f
#define ZOOM_SCALE_MIN 0.01f
#define ZOOM_SCALE_MAX 100.0f

#define LIGHTING_MIN 0.0f
#define LIGHTING_MAX 1.0f

#define SPECULAR_EXPONENT_LOG_MIN -2.0f
#define SPECULAR_EXPONENT_LOG_MAX 2.0f
#define SPECULAR_EXPONENT_MIN 0.01f
#define SPECULAR_EXPONENT_MAX 100.0f

// this is to scale to/from integer values
#define PARAMETER_SCALING 100


// end of include guard
#endif
