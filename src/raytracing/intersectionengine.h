#pragma once
#include <QList>
#include <raytracing/ray.h>
#include <scene/scene.h>
#include <raytracing/BVHnode.h>


class Scene;
class Ray;


class IntersectionEngine
{
public:
    IntersectionEngine();
    Scene *scene;

    BVHnode* BVHrootNode; //pointer to root of bvh tree

    Intersection GetIntersection(Ray r);
    QList<Intersection> GetAllIntersections(Ray r);
};
