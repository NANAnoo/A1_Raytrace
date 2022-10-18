#include "Triangle.h"
#include <math.h>
#include <limits>

Triangle::Triangle()
{
    valid = false;
    shared_material = nullptr;
}

float Triangle::intersect(Ray r)
{
    // TODO: Calculate intersection between a ray and a triangle!
    //  1. check if the ray is parallel to the triangle

    // get the normal of self, n = unit(V_1V_0 x V_2V_1)
    Cartesian3 planar_normal = ((verts[1] - verts[0]).Vector().cross((verts[2] - verts[1]).Vector())).unit();
    // if the ray.direction is perpendicular to the normal, then No intersect, return infinity
    if (std::abs(r.direction.dot(planar_normal)) < std::numeric_limits<float>::epsilon())
        return std::numeric_limits<float>::infinity();

    // 2. find the intersection P of the ray and the planar
    // For point A on planr, assume intersection is P
    // vector(A,P) * planar_normal = 0;
    // (vector(origin,P) - vector(origin,A)) * planar_normal = 0;
    // (direction * t - vector(origin,A)) * planar_normal = 0;
    // t = vector(origin,A) * planar_normal / (direction * planar_normal);
    float t = (verts[0] - r.origin).Vector().dot(planar_normal) / (r.direction.dot(planar_normal));
    if (t < 0)
        // the intersection point is on the opposite direction with the ray.
        return std::numeric_limits<float>::infinity();

    // 3. check if the point is in the triangle
    // check if the P is on the left side of each edge
    // check if the cross result has same direction with normal, if yes, the P is on left
    Cartesian3 intersection = r.origin + t * r.direction;
    if (((verts[1] - verts[0]).Vector().cross(intersection - verts[0].Vector())).dot(planar_normal) > 0 &&
        ((verts[2] - verts[1]).Vector().cross(intersection - verts[1].Vector())).dot(planar_normal) > 0 &&
        ((verts[0] - verts[2]).Vector().cross(intersection - verts[2].Vector())).dot(planar_normal) > 0)
        // P = origin + t * direction, use t to represent the intersection point
        return t;

    return std::numeric_limits<float>::infinity(); // just to compile warning free :)
}

void Triangle::validate()
{
    valid = true;
}

bool Triangle::isValid()
{
    return valid;
}

Cartesian3 Triangle::baricentric(Cartesian3 o)
{
    // TODO: Implement this! Input is the intersection between the ray and the triangle.

    // use the Si / S_tri to represent baricentric interpolation
    // S = ||vecA cross vecB|| / 2
    Cartesian3 e_01 = (verts[1] - verts[0]).Vector();
    Cartesian3 e_12 = (verts[2] - verts[1]).Vector();
    Cartesian3 e_20 = (verts[0] - verts[2]).Vector();
    Cartesian3 l_P0 = o - verts[0].Vector();
    Cartesian3 l_P1 = o - verts[1].Vector();
    // area of self;
    float area = (e_01.cross(e_12)).length();
    // area of three parts
    float area_of_triP12_times_2 = (l_P1.cross(e_12)).length();
    float area_of_tri0P2_times_2 = (l_P0.cross(e_20)).length();
    float area_of_tri01P_times_2 = (l_P0.cross(e_01)).length();
    
    Cartesian3 bc(area_of_triP12_times_2 / area, 
                  area_of_tri0P2_times_2 / area, 
                  area_of_tri01P_times_2 / area);
    if (bc.x + bc.y + bc.z < 0.99f) {
        std::cout<< "Oops, wrong !" << std::endl;
    }
    return bc;
}
