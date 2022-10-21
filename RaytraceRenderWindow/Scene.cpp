#include "Scene.h"
#include <limits>
#include <math.h>

Scene::Scene(std::vector<ThreeDModel> *texobjs, RenderParameters *renderp)
{
    objects = texobjs;
    rp = renderp;

    Cartesian3 ambient = Cartesian3(0.5f, 0.5f, 0.5f);
    Cartesian3 diffuse = Cartesian3(0.5f, 0.5f, 0.5f);
    Cartesian3 specular = Cartesian3(0.5f, 0.5f, 0.5f);
    Cartesian3 emissive = Cartesian3(0, 0, 0);
    float shininess = 1.0f;

    default_mat = new Material(ambient, diffuse, specular, emissive, shininess);
}

Matrix4 Scene::getModelview()
{
    Matrix4 result;
    // TODO: Grab all the necessary matrices to build your modelview. Sliders, Arcball, centering.
    Matrix4 transilation;
    transilation.SetTranslation(Cartesian3(rp->xTranslate, rp->yTranslate, rp->zTranslate - 1));
    // rotate and translate
    result = transilation * rp->rotationMatrix;
    return result;
}

// updateScene will build the scene like we would do for
// rasterization. Basically go through the model and
// create triangles. We however need to transform things
// to VCS to raytrace

void Scene::updateScene()
{
    triangles.clear();
    scene_lights.clear();

    Matrix4 modelView = getModelview();
    // We go through all the objects to construct the scene
    for (int i = 0; i < int(objects->size()); i++)
    {
        typedef unsigned int uint;
        ThreeDModel obj = objects->at(uint(i));

        // loop through the faces: note that they may not be triangles, which complicates life
        for (unsigned int face = 0; face < obj.faceVertices.size(); face++)
        { // per face
            // on each face, treat it as a triangle fan starting with the first vertex on the face
            for (unsigned int triangle = 0; triangle < obj.faceVertices[face].size() - 2; triangle++)
            { // per triangle
                // now do a loop over three vertices
                Triangle t;
                for (unsigned int vertex = 0; vertex < 3; vertex++)
                { // per vertex
                    // we always use the face's vertex 0
                    uint faceVertex = 0;
                    // so if it isn't 0, we want to add the triangle base ID
                    if (vertex != 0)
                        faceVertex = triangle + vertex;

                    // this is our vertex before any transformations. (world space)
                    Homogeneous4 v = Homogeneous4(
                        obj.vertices[obj.faceVertices[face][faceVertex]].x,
                        obj.vertices[obj.faceVertices[face][faceVertex]].y,
                        obj.vertices[obj.faceVertices[face][faceVertex]].z);

                    // This will start working when you write the method above!
                    v = modelView * v;

                    t.verts[vertex] = v;

                    Homogeneous4 n = Homogeneous4(
                        obj.normals[obj.faceNormals[face][faceVertex]].x,
                        obj.normals[obj.faceNormals[face][faceVertex]].y,
                        obj.normals[obj.faceNormals[face][faceVertex]].z,
                        0.0f);

                    n = modelView * n;
                    t.normals[vertex] = n;

                    Cartesian3 tex = Cartesian3(
                        obj.textureCoords[obj.faceTexCoords[face][faceVertex]].x,
                        obj.textureCoords[obj.faceTexCoords[face][faceVertex]].y,
                        0.0f);
                    t.uvs[vertex] = tex;

                    t.colors[vertex] = Cartesian3(0.7f, 0.7f, 0.7f);

                } // per vertex
                t.validate();
                if (obj.material == nullptr)
                {
                    t.shared_material = default_mat;
                }
                else
                {
                    t.shared_material = obj.material;
                }
                triangles.push_back(t);
            } // per triangle
        }     // per face
    }         // per object

    // transform all the lights in rp
    for (Light *l : rp->lights)
    {
        scene_lights.push_back(l->TransformedLight(modelView));
    }
}

Scene::CollisionInfo Scene::closestTriangle(Ray r)
{
    // TODO: method to find the closest triangle!

    Scene::CollisionInfo ci;
    ci.t = std::numeric_limits<float>::infinity();
    for (auto tri : triangles)
    {
        float cur_t = tri.intersect(r);
        // smaller t, closer triangle.
        if (ci.t > cur_t)
        {
            ci.t = cur_t;
            ci.tri = tri;
        }
    }

    return ci;
}

Homogeneous4 Scene::colorFromRay(Ray &r)
{
    Homogeneous4 color;

    std::vector<RecursiveInfo *> callStack;
    std::vector<RecursiveInfo *> returnStask;
    
    Scene::RecursiveInfo *root_info = new RecursiveInfo(r);
    root_info->return_compute_type = Scene::DoneType;
    root_info->material_from = nullptr;
    callStack.push_back(root_info);
    // recusive path_trace, shader frienly, no recursion function
    while(callStack.size() != 0) {
        RecursiveInfo *info = callStack.back();
        callStack.pop_back();
        if (returnStask.size() > 10) {
            // reach max depth
            info->need_compute_types = DoneType;
        } else {
            GetColorFromRayOneStep(info, callStack);
        }
        if (info->need_compute_types != Scene::DoneType) {
            returnStask.push_back(info);
        } else { 
            // trace back to compute father nodes color
            RecursiveInfo *returnInfo = nullptr;
            do {
                if (info == root_info) { 
                    // reach root result;
                    color = info->light_color.modulate(info->point_color);
                    delete info;
                    info == nullptr;
                    return color;
                }
                returnInfo = returnStask.back();
                if (returnInfo == nullptr) 
                    // somthing went wrong
                    return Homogeneous4(0, 0, 0, 0);
                float returnMaterial;
                switch (info->return_compute_type)
                {
                    case ReflecType:
                        returnMaterial = info->material_from->reflectivity;
                        break;
                    case RefractType:
                        returnMaterial = info->material_from->transparency;   
                    default:
                        break;
                }
                returnInfo->light_color = returnInfo->light_color + info->light_color.modulate(info->point_color) * returnMaterial;
                returnInfo->need_compute_types = static_cast<ComputeType>(returnInfo->need_compute_types - info->return_compute_type);
                delete info;
                info = nullptr;
                if (returnInfo->need_compute_types == DoneType) { 
                    // continue trace back
                    returnStask.pop_back();
                    info = returnInfo;
                }
            } while (info != nullptr);
        }
    }
    
    return color;
}

void Scene::GetColorFromRayOneStep(Scene::RecursiveInfo *info, std::vector<Scene::RecursiveInfo *> &callstack)
{
    Homogeneous4 color;
    Ray r = info->ray_from;
    Scene::CollisionInfo ci = closestTriangle(r);
    Cartesian3 intersetion = r.origin + ci.t * r.direction;
    // calculate color provided from the light;
    if (ci.t <= std::numeric_limits<float>::epsilon() || ci.t >= std::numeric_limits<float>::infinity())
        return;

    // get interpolation normal and color
    Cartesian3 baricentric = ci.tri.baricentric(intersetion);
    Cartesian3 normal = (baricentric.x * ci.tri.normals[0].Vector() +
                        baricentric.y * ci.tri.normals[1].Vector() +
                        baricentric.z * ci.tri.normals[2].Vector()).unit();

    Homogeneous4 tri_color = baricentric.x * ci.tri.colors[0] +
                             baricentric.y * ci.tri.colors[1] +
                             baricentric.z * ci.tri.colors[2];
    info->point_color = tri_color;
    if (ci.tri.shared_material->isLight())
    {
        info->light_color = ci.tri.shared_material->emissive;
        info->point_color = Homogeneous4(1, 1, 1, 1);
        info->need_compute_types = DoneType;
        return;
    }
    intersetion = intersetion + std::numeric_limits<float>::epsilon() * normal;
    // Blihn-Phong & Shade from lights:
    if ((rp->phongEnabled || rp->shadowsEnabled) && (r.ray_type == Ray::primary))
    {
        for (Light *l : scene_lights)
        {
            if (rp->shadowsEnabled && CheckShadowATPoint(intersetion, normal, l))
            {
                // light is blocked, if there is Blihn-Phong, add ambient
                if (rp->phongEnabled)
                    color = color + tri_color.modulate(ci.tri.shared_material->ambient);
                continue;
            }
            // Blihn-Phong
            if (rp->phongEnabled)
                color = color + tri_color.modulate(GetColorFromBlinnPhongAtPoint(r.origin, intersetion, normal, ci.tri.shared_material, l));
            else
                // default color
                color = color + tri_color.modulate(ci.tri.shared_material->ambient);
        }
    }
    if (rp->reflectionEnabled && ci.tri.shared_material->reflectivity > 0.f) {
        info->need_compute_types = static_cast<ComputeType>(info->need_compute_types + ReflecType);
        RecursiveInfo *newInfo = new RecursiveInfo(r.getReflectAt(intersetion, normal));
        newInfo->return_compute_type = ReflecType;
        newInfo->material_from = ci.tri.shared_material;
        callstack.push_back(newInfo);
    }
    if (rp->refractionEnabled && ci.tri.shared_material->transparency > 0.f) {
        bool valid;
        Ray refract_ray = r.getRefractAt(intersetion, normal, info->previous_ior, ci.tri.shared_material->indexOfRefraction, valid);
        if (valid) {
            info->need_compute_types = static_cast<ComputeType>(info->need_compute_types + RefractType);
            RecursiveInfo *newInfo = new RecursiveInfo(refract_ray);
            newInfo->return_compute_type = RefractType;
            newInfo->previous_ior = info->material_from == nullptr ? 1.f : info->material_from->indexOfRefraction;
            newInfo->material_from = ci.tri.shared_material;
            callstack.push_back(newInfo);
        }
        if (rp->fresnelRendering) {
            // fresnel reflect

        }
    }
    if (rp->monteCarloEnabled) {

    }
    info->light_color = color;
}

Homogeneous4 Scene::GetColorFromBlinnPhongAtPoint(Cartesian3 &lookFrom, Cartesian3 &hitPoint, Cartesian3 &normal, Material *m, Light *l)
{
    Cartesian3 light_position = l->GetPositionCenter().Vector();
    Cartesian3 O_L = light_position - hitPoint;
    Cartesian3 U_O_L = O_L.unit();
    Cartesian3 U_O_R = (lookFrom - hitPoint).unit();
    float length_OL = O_L.length();
    // half_angle_between_normal_and_centerline
    float angle_half = std::max(0.f, (U_O_L + U_O_R).unit().dot(normal));
    // angle_of_normal_and_light
    float angle_nl = std::max(0.f, normal.dot(U_O_L));
    // distance weaken
    float distance_factor = 1.f / (length_OL * length_OL);

    Homogeneous4 specular_light = m->specular * powf(angle_half, m->shininess);
    Homogeneous4 diffuse_light = m->diffuse * angle_nl;
    // add ambient, specular_light, diffuse_light together
    return Homogeneous4(m->ambient) + distance_factor * (specular_light + diffuse_light).modulate(l->GetColor());
}

bool Scene::CheckShadowATPoint(Cartesian3 &hitPoint, Cartesian3 &normal, Light *l)
{
    Cartesian3 lightCenter = l->GetPositionCenter().Vector();
    Ray lightr(hitPoint, (lightCenter - hitPoint).unit(), Ray::Type::primary);
    float EPS = std::numeric_limits<float>::epsilon();
    // To avoid the light hitting the intersetion's triangle
    // float dis = 10 * EPS / lightr.direction.dot(normal);
    // lightr.origin = lightr.origin + lightr.direction * dis;

    CollisionInfo tci = closestTriangle(lightr);
    float t_max = (lightCenter - lightr.origin).length();
    return tci.t > EPS && tci.t < t_max && !tci.tri.shared_material->isLight();
}
