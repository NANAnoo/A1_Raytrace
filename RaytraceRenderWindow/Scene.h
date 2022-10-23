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
    std::vector<Light> scene_lights;
    Scene(std::vector<ThreeDModel> *texobjs, RenderParameters *renderp);
    void updateScene();
    CollisionInfo closestTriangle(Ray r);
    Matrix4 getModelview();

    // main render funtion
    // Homogeneous4 colorFromRay(Ray &r);
    Homogeneous4 colorFromRay(Ray &r, float current_IOR, int max_depth);

private:
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
    Homogeneous4 GetColorFromBlinnPhongAtPoint(Cartesian3 &lookFrom, Cartesian3 &hitPoint, Cartesian3 &normal, Material *m, Light l);

    /**
     * @brief check if there is a shadow at this point from a light
     *
     * @param hitPoint
     * @param normal
     * @param l
     * @return true
     * @return false
     */
    bool CheckShadowATPoint(Cartesian3 &hitPoint, Cartesian3 &light_position);

    // optimize, AABB box & BVH 
    struct AABB
    {
        Cartesian3 min;
        Cartesian3 max;
        Triangle tri;
        bool hit(const Ray &r, float tmin, float tmax);
        AABB() : min(Cartesian3()), max(Cartesian3()){}
        AABB(const AABB & other) : min(other.min), max(other.max), tri(other.tri){}
        AABB(Triangle &tri);
        AABB(const Cartesian3 &l, const Cartesian3 &r) : min(l), max(r){}
        const AABB operator=(const AABB other);
    };

    struct BVHNode
    {
        BVHNode *left;
        BVHNode *right;
        AABB boundingBox;
        BVHNode(AABB &box) : boundingBox(box), left(nullptr), right(nullptr){}
        BVHNode(std::vector<AABB> &boxes, unsigned int begin, unsigned int end);
        bool hit(Ray r, float tmin, float tmax, CollisionInfo &ci);
    };
    
    BVHNode *bvh_root;
};

#endif // SCENE_H
