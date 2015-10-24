#include <scene/geometry/geometry.h>

float Geometry::RayPDF(const Intersection &isx, const Ray &ray)
{
    //The isx passed in was tested ONLY against us (no other scene objects), so we test if NULL
    //rather than if != this.
    if(isx.object_hit == NULL)
    {
        return 0;
    }

//find THETA between incoming_ray and surface normal at pt.

    //normalize both the ray and the normal
    glm::vec3 incoming_ray = glm::normalize(ray.direction);
    glm::vec3 norm = glm::normalize(isx.normal);
    //take dot product of these vectors
    float val = glm::dot(incoming_ray, norm);
    //solve for THETA
    float theta = glm::acos(val);

    return 1.0f/(glm::cos(theta)*this->area);

}

Intersection Geometry::SamplePoint(float a, float b) {
    Intersection isx;
    return isx;
}
