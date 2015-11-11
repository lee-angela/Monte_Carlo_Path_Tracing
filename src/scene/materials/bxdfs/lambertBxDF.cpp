#include <scene/materials/bxdfs/lambertBxDF.h>

glm::vec3 LambertBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //integrating over entire hemisphere of visible directions from p
    //integral will eval to base_color
    return this->diffuse_color/PI;
}


glm::vec3 BxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    //sample random point on hemisphere using Malley's method:
    glm::vec4 pt1 = glm::vec4(0.5f,0.0f,0.0f,1.0f); //any point on circumfrence
    float angle = rand1*2*PI; //fraction of the entire circumfrence of disc
    glm::vec3 axis(0,0,1); //rotate around z axis

    glm::vec3 pt2 = glm::vec3(glm::rotate(glm::mat4(1.0f), angle, axis) * pt1);
    glm::vec4 discPt = glm::vec4(rand2*pt2,1.0f); //between pt2 and origin in Cartesian coordinates IN LOCAL SPACE (origin is at 0,0,0)
    //convert to spherical coords
    float p = glm::sqrt(glm::pow(discPt[0],2)+glm::pow(discPt[1],2)+glm::pow(discPt[2],2));
    float phi = glm::acos(discPt[2]/p);
    glm::vec3 sph_pt = glm::vec3(p, phi, 0.0f); //spherical coords of random point on disc

    //project onto hemisphere
//    float theta = glm::asin(p);

    //calculate hemisphere sampled pt:
    discPt[0] = p*glm::cos(phi);
    discPt[1] = p*glm::sin(phi);
    discPt[2] = glm::sqrt(1-(glm::pow(discPt[0],2) + glm::pow(discPt[1],2)));

    //find hemisphere sampled pt & ray, send to EvaluateScatteredEnergy()
    glm::vec3 energy = EvaluateScatteredEnergy(wo,discPt);

    *wi_ret = discPt;
    pdf_ret = 1/PI;
    return energy;
}
