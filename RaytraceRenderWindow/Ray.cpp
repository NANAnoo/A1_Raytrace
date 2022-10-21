#include "Ray.h"

Ray::Ray(Cartesian3 og, Cartesian3 dir, Type rayType)
{
    origin = og;
    direction = dir;
    ray_type = rayType;
}

Ray Ray::getReflectAt(Cartesian3 &hit, Cartesian3 &normal)
{
    Cartesian3 retDirection = direction + normal * 2;
    return Ray(hit, retDirection, secondary);
}

Ray Ray::getRefractAt(Cartesian3 &hit, Cartesian3 normal, float t)
{

}

Ray Ray::getRandomReflect(Cartesian3 &hit)
{

}
