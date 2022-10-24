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
    bvh_root = nullptr;

    default_mat = new Material(ambient, diffuse, specular, emissive, shininess);
}

Scene::~Scene()
{
    destoryBVHTree(bvh_root);
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
    
    // destory tree if need
    destoryBVHTree(bvh_root);
    // build BVH tree
    if (triangles.size() > 10)
    {
        std::vector<AABB> boxes;
        for (Triangle &tri : triangles)
        {
            boxes.push_back(tri);
        }
        // create bvh tree
        bvh_root = new BVHNode(boxes, 0, boxes.size() - 1);
    }
}

Scene::CollisionInfo Scene::closestTriangle(Ray r)
{
    // TODO: method to find the closest triangle!
    if (bvh_root != nullptr)
    {
        CollisionInfo ci;
        bvh_root->hit(r, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::infinity(), ci);
        return ci;
    }

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

Homogeneous4 Scene::colorFromRay(Ray &r, float current_IOR, int maxdepth)
{
    Homogeneous4 color(0, 0, 0, 0);
    if (maxdepth <= 0)
    {
        return color;
    }
    Scene::CollisionInfo ci = closestTriangle(r);
    // calculate color provided from the light;
    if (ci.t <= std::numeric_limits<float>::epsilon() || ci.t >= std::numeric_limits<float>::infinity() || !ci.tri.isValid())
        return color;

    // get interpolation normal and color
    Cartesian3 intersetion = r.origin + ci.t * r.direction;
    Cartesian3 baricentric = ci.tri.baricentric(intersetion);
    Cartesian3 normal = (baricentric.x * ci.tri.normals[0].Vector() +
                         baricentric.y * ci.tri.normals[1].Vector() +
                         baricentric.z * ci.tri.normals[2].Vector())
                            .unit();
    if (rp->interpolationRendering)
    {
        return (normal + Cartesian3(1, 1, 1)) / 2;
    }

    Homogeneous4 tri_color = baricentric.x * ci.tri.colors[0] +
                             baricentric.y * ci.tri.colors[1] +
                             baricentric.z * ci.tri.colors[2];
    Material *m = ci.tri.shared_material;

    if (m->texture != nullptr)
    {
        Cartesian3 texCoord = baricentric.x * ci.tri.uvs[0] +
                              baricentric.y * ci.tri.uvs[1] +
                              baricentric.z * ci.tri.uvs[2];
        RGBAValue c = m->texture->GetTexel(texCoord.x, texCoord.y, true);
        tri_color.w = c.alpha / 255.f;
        tri_color.x = c.red / 255.f;
        tri_color.y = c.green / 255.f;
        tri_color.z = c.blue / 255.f;
    }

    if (m->isLight())
    {
        return m->emissive;
    }
    float EPS = std::numeric_limits<float>::epsilon();
    intersetion = intersetion + 10 * EPS * normal;
    // indirectcolor
    // amblight
    bool needPhong = true;
    float cosine = r.direction.dot(normal);
    if (rp->monteCarloEnabled && (1.f - m->reflectivity - m->transparency) > EPS && cosine < 0)
    {
        // enable monteCarlo
        Ray difuzz_ray = r.getRandomReflect(intersetion, normal);
        difuzz_ray.ray_type = r.ray_type;
        float ambient = (1.f - m->transparency - m->reflectivity);
        Homogeneous4 temp_color = ambient * colorFromRay(difuzz_ray, current_IOR, maxdepth - 1).modulate(m->ambient);
        // pdf = 1 / (2 * pi)
        color = color + 2.f * static_cast<float>(M_PI) * temp_color;
        needPhong = false;
        //        for (Light l : scene_lights) {
        //            Cartesian3 light_position = l.GetPosition().Vector();
        //            if (CheckShadowATPoint(intersetion, light_position))
        //                continue;
        //            Cartesian3 O_L = light_position - intersetion;
        //            float r = O_L.length();
        //            float cosine = std::max(0.f, normal.dot(O_L));
        //            color = color + (1.f / (r * r) * cosine * l.GetColor()).modulate(m->diffuse);
        //        }
    }
    // reflect
    if (rp->reflectionEnabled && m->reflectivity > EPS && cosine < 0)
    {
        Ray reflect_ray = r.getReflectAt(intersetion, normal);
        reflect_ray.ray_type = r.ray_type;
        color = color + m->reflectivity * colorFromRay(reflect_ray, current_IOR, maxdepth - 1);
        needPhong = false;
    }

    // transparency
    float refract_part = 1.f;
    float ior_from = current_IOR, ior_to = m->indexOfRefraction;
    bool go_out = r.direction.dot(normal) > 0;
    if (go_out)
    {
        ior_from = ior_to;
        ior_to = 1.f;
    }
    // Fresnel
    if (rp->fresnelRendering && (m->transparency > EPS || !rp->monteCarloEnabled))
    {
        float prob = r.getFresnelProbability(normal, ior_from, ior_to);
        refract_part = 1.f - prob;
        Cartesian3 fresnel_normal = go_out ? -1.f * normal : normal;
        Ray fresnel_ray = r.getReflectAt(intersetion, fresnel_normal);
        fresnel_ray.ray_type = r.ray_type;
        if (r.ray_type == Ray::primary)
        {
            // for internal refraction
            fresnel_ray.origin = fresnel_ray.origin + 20.f * EPS * fresnel_normal;
        }
        if (m->transparency > EPS)
            color = color + m->transparency * prob * colorFromRay(fresnel_ray, current_IOR, maxdepth - 1);
        else
            color = color + prob * colorFromRay(fresnel_ray, current_IOR, maxdepth - 1).modulate(m->ambient);
    }
    // refraction
    if (rp->refractionEnabled && m->transparency > EPS)
    {
        bool valid;
        Ray refract_ray = r.getRefractAt(intersetion, normal, ior_from, ior_to, valid);
        if (valid)
        {
            color = color + m->transparency * refract_part * colorFromRay(refract_ray, ior_to, maxdepth - 1);
        }
        needPhong = false;
    }
    // Blihn-Phong & Shade from lights:
    // direct color
    if (needPhong && (rp->phongEnabled || rp->shadowsEnabled))
    {
        Homogeneous4 phong_color(0, 0, 0, 0);
        for (Light &l : scene_lights)
        {
            // disabled under monteCarlo
            Cartesian3 light_position = l.GetPositionCenter().Vector();
            if (rp->shadowsEnabled && CheckShadowATPoint(intersetion, light_position))
                continue;
            // Blihn-Phong
            if (rp->phongEnabled)
                phong_color = phong_color + GetColorFromBlinnPhongAtPoint(r.origin, intersetion, normal, m, l);
            else
                // default color
                color = Homogeneous4(1, 1, 1, 1);
        }
        if (rp->phongEnabled)
            phong_color = phong_color + m->ambient;
        color = color + tri_color.modulate(phong_color);
    }
    return color;
}

Homogeneous4 Scene::GetColorFromBlinnPhongAtPoint(Cartesian3 &lookFrom, Cartesian3 &hitPoint, Cartesian3 &normal, Material *m, Light l)
{
    Cartesian3 light_position = l.GetPositionCenter().Vector();
    Cartesian3 O_L = light_position - hitPoint;
    Cartesian3 U_O_L = O_L.unit();
    Cartesian3 U_O_R = (lookFrom - hitPoint).unit();
    float length_OL = O_L.length();
    // half_angle_between_normal_and_centerline
    float angle_half = std::max(0.f, (U_O_L + U_O_R).unit().dot(normal));
    // angle_of_normal_and_light
    float angle_nl = std::max(0.f, normal.dot(U_O_L));
    // distance weaken
    float distance_factor = std::min(1.f, 1.f / (length_OL * length_OL));

    Homogeneous4 specular_light = m->specular * powf(angle_half, m->shininess);
    Homogeneous4 diffuse_light = m->diffuse * angle_nl;
    // add ambient, specular_light, diffuse_light together
    return distance_factor * (specular_light + diffuse_light).modulate(l.GetColor());
}

bool Scene::CheckShadowATPoint(Cartesian3 &hitPoint, Cartesian3 &light_position)
{
    Ray lightr(hitPoint, (light_position - hitPoint).unit(), Ray::Type::primary);
    float EPS = std::numeric_limits<float>::epsilon();
    // To avoid the light hitting the intersetion's triangle

    CollisionInfo tci = closestTriangle(lightr);
    float t_max = (light_position - lightr.origin).length();
    return tci.t > EPS && tci.t < t_max && !tci.tri.shared_material->isLight() && !(tci.tri.shared_material->transparency > 0.f && rp->refractionEnabled);
}

Scene::AABB::AABB(Triangle &tri)
{
    float EPS = std::numeric_limits<float>::epsilon();
    min = Cartesian3(
        fmin(fmin(tri.verts[0][0], tri.verts[1][0]), tri.verts[2][0]) - EPS,
        fmin(fmin(tri.verts[0][1], tri.verts[1][1]), tri.verts[2][1]) - EPS,
        fmin(fmin(tri.verts[0][2], tri.verts[1][2]), tri.verts[2][2]) - EPS);
    max = Cartesian3(
        fmax(fmax(tri.verts[0][0], tri.verts[1][0]), tri.verts[2][0]) + EPS,
        fmax(fmax(tri.verts[0][1], tri.verts[1][1]), tri.verts[2][1]) + EPS,
        fmax(fmax(tri.verts[0][2], tri.verts[1][2]), tri.verts[2][2]) + EPS);
    this->tri = tri;
}

const Scene::AABB Scene::AABB::operator=(const Scene::AABB other)
{
    this->min = other.min;
    this->max = other.max;
    this->tri = other.tri;
    return *this;
}

bool Scene::AABB::hit(const Ray &r, float tmin, float tmax)
{
    for (int i = 0; i < 3; i++)
    {
        float t0 = fmin(
            (min[i] - r.origin[i]) / r.direction[i],
            (max[i] - r.origin[i]) / r.direction[i]);
        float t1 = fmax(
            (min[i] - r.origin[i]) / r.direction[i],
            (max[i] - r.origin[i]) / r.direction[i]);
        tmin = fmax(t0, tmin);
        tmax = fmin(t1, tmax);
        if (tmax < tmin)
        {
            return false;
        }
    }
    return true;
}

Scene::BVHNode::BVHNode(std::vector<AABB> &boxes, unsigned int begin, unsigned int end)
{
    // build tree
    int size = end - begin + 1;
    switch (size)
    {
    case 1:
        left = new BVHNode(boxes[begin]);
        right = nullptr;
        break;
    case 2:
        left = new BVHNode(boxes[begin]);
        right = new BVHNode(boxes[end]);
        break;
    default:
        // sort objs
        AABB curmin;
        for (unsigned int i = begin; i <= end; i++)
        {
            curmin = boxes[i];
            unsigned int index = i;
            for (unsigned int j = i + 1; j <= end; j++)
            {
                AABB other = boxes[j];
                int r = rand() % 3;
                if (curmin.min[r] > other.min[r])
                {
                    curmin = other;
                    index = j;
                }
            }
            if (index != i) {
                boxes[index] = boxes[i];
                boxes[i] = curmin;
            }
        }
        unsigned int mid = (begin + end) / 2;
        left = new BVHNode(boxes, begin, mid);
        right = new BVHNode(boxes, mid + 1, end);
        break;
    }
    AABB l = left->boundingBox, r;
    if (right != nullptr)
        r = right->boundingBox;
    else
        r = l;
    // union two boxes
    boundingBox.min = Cartesian3(
        fmin(l.min.x, r.min.x),
        fmin(l.min.y, r.min.y),
        fmin(l.min.z, r.min.z));
    boundingBox.max  = Cartesian3(
        fmax(l.max.x, r.max.x),
        fmax(l.max.y, r.max.y),
        fmax(l.max.z, r.max.z));
}

bool Scene::BVHNode::hit(Ray r, float tmin, float tmax, CollisionInfo &ci)
{

    if (boundingBox.hit(r, tmin, tmax))
    {
        if (left == nullptr && right == nullptr)
        {
            if (boundingBox.tri.isValid())
            {
                float t = boundingBox.tri.intersect(r);
                ci.tri = boundingBox.tri;
                ci.t = t;
                if (t > tmin && t < tmax)
                {
                    return true;
                }
            }
            ci.t = std::numeric_limits<float>::infinity();
            return false;
        }
        CollisionInfo l_ci, r_ci;
        bool hit_left = left->hit(r, tmin, tmax, l_ci);
        bool hit_right = false;
        if (right != nullptr)
            hit_right = right->hit(r, tmin, tmax, r_ci);
        if (hit_left && hit_right)
        {
            ci = l_ci.t < r_ci.t ? l_ci : r_ci;
            return true;
        }
        if (hit_left)
        {
            ci = l_ci;
            return true;
        }
        if (hit_right)
        {
            ci = r_ci;
            return true;
        }
    }
    ci.t = std::numeric_limits<float>::infinity();
    return false;
}


void Scene::destoryBVHTree(BVHNode *node)
{
    if (node == nullptr) {
        return;
    }
    destoryBVHTree(node->left);
    destoryBVHTree(node->right);
    delete node;
}