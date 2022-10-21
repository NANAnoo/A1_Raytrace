#include "Ray.h"

#include <math.h>

Ray::Ray(Cartesian3 og, Cartesian3 dir, Type rayType)
{
    origin = og;
    direction = dir;
    ray_type = rayType;
}

Ray Ray::getReflectAt(Cartesian3 &hit, Cartesian3 &normal)
{
    Cartesian3 retDirection = direction + normal * 2;
    return Ray(hit, retDirection, primary);
}

Ray Ray::getRefractAt(Cartesian3 &hit, Cartesian3 &normal, float ior_from, float ior_to, bool &valid)
{
    // sin(theta_to) = sin(theta_from) * n1 / n2
    float cos_theta_from, sin_theta_from, cos_theta_to, sin_theta_to;
    bool swaped = false;
    if ((cos_theta_from = normal.dot(direction)) > 0) {
        // go out
        std::swap(ior_from, ior_to);
        swaped = true;
    } else {
        cos_theta_from = - cos_theta_from;
    }
    sin_theta_from = sqrtf(1.f - cos_theta_from * cos_theta_from);
    sin_theta_to = sin_theta_from * ior_from / ior_to;
    if (sin_theta_to > 1) {
        valid = false;
        return Ray(Cartesian3(), Cartesian3(), primary);
    }
    cos_theta_to = sqrtf(1.f - sin_theta_to * sin_theta_to);
    // local coord system
    valid = true;
    Cartesian3 i = (swaped ? normal * cos_theta_from : -1.f * normal * cos_theta_from);
    Cartesian3 j = (direction - i).unit();
    Cartesian3 returnDirection = i.unit() * cos_theta_to + j * sin_theta_to;
    //float lift_step = swaped ? std::numeric_limits<float>::epsilon() : - std::numeric_limits<float>::epsilon();
    return Ray(hit, returnDirection, swaped ? primary : secondary);
}

Ray Ray::getRandomReflect(Cartesian3 &hit)
{
    return Ray(Cartesian3(), Cartesian3(), primary);
}
