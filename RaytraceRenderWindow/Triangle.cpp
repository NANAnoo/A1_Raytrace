#include "Triangle.h"
#include <math.h>

Triangle::Triangle()
{
    valid = false;
    shared_material= nullptr;
}

float Triangle::intersect(Ray r)
{
    //TODO: Calculate intersection between a ray and a triangle!
    return r.origin.x*0; // just to compile warning free :)
}

void Triangle::validate(){
    valid = true;
}

bool Triangle::isValid(){
    return valid;
}

Cartesian3 Triangle::baricentric(Cartesian3 o)
{
   //TODO: Implement this! Input is the intersection between the ray and the triangle.
   Cartesian3 bc;
   bc.x = o.x ; // Just to compile warning free :)
   return bc;
}
