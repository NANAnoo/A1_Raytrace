#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Homogeneous4.h"
#include "Ray.h"
#include "Material.h"
#include "RGBAImage.h"

class Triangle
{
private:
    bool valid;

public:
    Homogeneous4 verts[3];
    Homogeneous4 normals[3];
    Homogeneous4 colors[3];
    Cartesian3 uvs[3];

    Material *shared_material;

    Triangle();
    void validate();
    bool isValid();
    float intersect(Ray r);

    Cartesian3 baricentric(Cartesian3 o);


};

#endif // TRIANGLE_H
