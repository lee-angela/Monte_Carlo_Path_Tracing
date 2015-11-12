#include <scene/materials/bxdfs/lambertBxDF.h>

glm::vec3 LambertBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //integrating over entire hemisphere of visible directions from p
    //integral will eval to base_color
    return this->diffuse_color/PI;
}


glm::vec3 BxDF::SampleAndEvaluateScatteredEnergy(Intersection isx, const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    //sample random point on hemisphere using Malley's method:
    glm::vec3 discPt = UnifSampleDisc(rand1, rand2);
    //convert to spherical coords
    float p = glm::sqrt(glm::pow(discPt[0],2)+glm::pow(discPt[1],2)+glm::pow(discPt[2],2));
    float phi = glm::acos(discPt[2]/p);
    glm::vec3 sph_pt = glm::vec3(p, phi, 0.0f); //spherical coords of random point on disc
    //project onto hemisphere
    float theta = glm::asin(p);

    //calculate hemisphere sampled pt in Cartesian coords:
    discPt[0] = p*glm::cos(phi);
    discPt[1] = p*glm::sin(phi);
    discPt[2] = glm::sqrt(1-(glm::pow(discPt[0],2) + glm::pow(discPt[1],2)));

    //find hemisphere sampled pt & ray, send to EvaluateScatteredEnergy()
    glm::vec3 energy = EvaluateScatteredEnergy(wo,discPt);

    wi_ret = discPt;
    pdf_ret = glm::cos(theta)/PI;
    return energy;
}

glm::vec3 UnifSampleDisc(float u1, float u2)
{
    float sx = 2 * u1 - 1.0f;
    float sy = 2 * u2 - 1.0f;
    float r, theta;

    if (sx == 0.0 && sy == 0.0)
    {
        return glm::vec3(0,0,0);
    }
    if (sx >= -sy)
    {
        if (sx > sy)
        {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0f + sy/r;
        }
        else
        {
            // Handle second region of disk
            r = sy;
            theta = 2.0f - sx/r;
        }
    }
    else
    {
        if (sx <= sy)
        {
            // Handle third region of disk
            r = -sx;
            theta = 4.0f - sy/r;
        }
        else
        {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= PI / 4.f;
    float dx = r * cosf(theta);
    float dy = r * sinf(theta);
    return glm::vec3(dx, dy, 0);
}
