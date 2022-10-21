//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2022
//
//  ------------------------
//  Scene.h
//  ------------------------
//
//  Contains a definition of a scene, with triangles and transformations.
//
///////////////////////////////////////////////////

#ifndef SCENE_H
#define SCENE_H

#include "Homogeneous4.h"
#include "ThreeDModel.h"
#include <vector>
#include "Ray.h"
#include "Triangle.h"
#include "Material.h"

class Scene
{
public:
    struct CollisionInfo
    {
        Triangle tri;
        float t;
    };

    std::vector<ThreeDModel> *objects;
    RenderParameters *rp;
    Material *default_mat;

    std::vector<Triangle> triangles;
    std::vector<Light *> scene_lights;
    Scene(std::vector<ThreeDModel> *texobjs, RenderParameters *renderp);
    void updateScene();
    CollisionInfo closestTriangle(Ray r);
    Matrix4 getModelview();

    // main render funtion
    Homogeneous4 colorFromRay(Ray &r);

private:
    enum ComputeType
    {
        DoneType = 0 << 0,
        ReflecType = 1 << 0,
        RefractType = 1 << 8,
        AmbientType = 1 << 16, // may have more than one, put into last
    };
    struct RecursiveInfo
    {
        /* data */
        Ray ray_from;
        Material *material_from;
        Homogeneous4 point_color;
        Homogeneous4 light_color;
        ComputeType need_compute_types;
        ComputeType return_compute_type;
        float previous_ior;
        explicit RecursiveInfo(Ray r) :ray_from(r), need_compute_types(DoneType), return_compute_type(DoneType), previous_ior(1.f){}
    };
    /**
     * @brief Get the Color From Blinn Phong At Hit Point
     *
     * @param lookFrom
     * @param hitPoint
     * @param normal
     * @param m
     * @param l
     * @return light color
     */
    Homogeneous4 GetColorFromBlinnPhongAtPoint(Cartesian3 &lookFrom, Cartesian3 &hitPoint, Cartesian3 &normal, Material *m, Light *l);

    /**
     * @brief check if there is a shadow at this point from a light
     *
     * @param hitPoint
     * @param normal
     * @param l
     * @return true
     * @return false
     */
    bool CheckShadowATPoint(Cartesian3 &hitPoint, Cartesian3 &normal, Light *l);

    /**
     * @brief Get the Color From Ray in One Step
     *
     * @param info recusive ray info
     * @param callstack recusive stack
     */
    void GetColorFromRayOneStep(RecursiveInfo *info, std::vector<RecursiveInfo *> &callstack);

    /**
     * @brief Get the Reflect Ray from input
     *
     * @param input input ray
     * @param normal normal at hit point
     * @return Ray
     */
    Ray GetReflectRay(Ray &input, Cartesian3 &normal);

    /**
     * @brief Get the Refract Ray object
     *
     * @param input input ray
     * @param normal normal at hit point
     * @param n factor
     * @return Ray
     */
    Ray GetRefractRay(Ray &input, Cartesian3 &normal, float n);

    Ray GetAmbientRay(Ray &input, Cartesian3 &normal);
};

#endif // SCENE_H
