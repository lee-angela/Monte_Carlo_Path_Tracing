#pragma once

#include <scene/geometry/geometry.h>

//A cube is assumed to have side lengths of 1 and a center of <0,0,0>. This means all vertices are of the form <+/-0.5, +/-0.5, +/-0.5>
//These attributes can be altered by applying a transformation matrix to the cube.
class Cube : public Geometry
{
public:
    Intersection GetIntersection(Ray r);
    virtual glm::vec2 GetUVCoordinates(const glm::vec3 &point);
    virtual glm::vec3 ComputeNormal(const glm::vec3 &P);
    void create();
    void setBoundingBox();

    virtual void ComputeArea();

    //not implemented as light source
    //isx = Intersection on non-light geom
    virtual Intersection SamplePoint(float a, float b, Intersection isx);
    bool isMesh();
};
