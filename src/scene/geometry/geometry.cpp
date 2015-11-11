#include <scene/geometry/geometry.h>

/**
 * @brief Geometry::RayPDF
 * @param isx with light
 * @param ray going towards THIS geom
 * @return
 */
float Geometry::RayPDF(const Intersection &isx, const Ray &ray)

{
    //The isx passed in was tested ONLY against us (no other scene objects), so we test if NULL
    //rather than if != this.
    if(isx.object_hit == NULL)
    {
        return 0;
    }

    //ray leaving Geometry (leaving light)
    //isx is on light

    //ray leaving the light source
    glm::vec3 wj = glm::normalize(ray.origin - isx.point); // towards geom (in same dir as light normal)
    glm::vec3 norm = glm::normalize(isx.normal);
    //take dot product of these vectors
    float val = glm::dot(wj, norm);
    //solve for THETA
    float theta = glm::acos(val);
    float r = glm::distance(isx.point,ray.origin); //dist between light point & geom intersected
    return r*r/(glm::cos(theta)*this->area);

}

Intersection Geometry::SamplePoint(float a, float b, Intersection isx) {
    Intersection pt;
    return pt;
}
