#include "Ray.h"

#include <math.h>
#include <limits>

const float MY_PI = static_cast<float>(std::acos(-1));

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
    bool go_out = false;
    if ((cos_theta_from = normal.dot(direction)) > 0) {
        go_out = true;
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
    Cartesian3 i = (go_out ? normal * cos_theta_from : -1.f * normal * cos_theta_from);
    Cartesian3 j = (direction - i).unit();
    i = i.unit();
    Cartesian3 returnDirection = i * cos_theta_to + j * sin_theta_to;
    if (go_out) {
         return Ray(hit, returnDirection, go_out ? primary : secondary);
    }
    return Ray(hit + 20.f * std::numeric_limits<float>::epsilon() * i, returnDirection, go_out ? primary : secondary);
}

float Ray::getFresnelProbability(Cartesian3 &normal, float ior_from, float ior_to)
{
    float ior_rate = ior_from / ior_to, cos_theta_from;
    if ((cos_theta_from = normal.dot(direction)) > 0) {
        // go out
        ior_rate = 1.f / ior_rate;
    } else {
        cos_theta_from = - cos_theta_from;
    }
    float r0 = (1.f - ior_rate) / (1.f + ior_rate);
    r0 = r0 * r0;
    return (r0 + (1 - r0) * powf((1 - cos_theta_from), 5));
}

Ray Ray::getRandomReflect(Cartesian3 &hit, Cartesian3 &normal)
{
    float alpha = (static_cast<float>(std::rand()) / RAND_MAX) * MY_PI;
    float theta = (static_cast<float>(std::rand()) / RAND_MAX) * MY_PI * 2;
    Cartesian3 direct(std::sin(alpha) * std::cos(theta),std::cos(alpha),std::sin(alpha) * std::sin(theta));
    float cosine = direct.dot(normal);
    if (cosine < 0) {
        direct = direct - 2 * cosine * normal;
    }
    return Ray(hit, direct, primary);
}
