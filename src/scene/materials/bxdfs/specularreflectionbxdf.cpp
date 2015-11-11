#include <scene/materials/bxdfs/specularreflectionBxDF.h>

glm::vec3 SpecularReflectionBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //should be light color?
    return this->reflection_color; //perfectly specular will always reflect in direction of sampled wi!
}


glm::vec3 BxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    glm::vec3 norm = glm::vec3(0.0f,1.0f,0.0f);
    glm::vec3 wo_norm = glm::normalize(wo);
    wi_ret = glm::normalize(wo_norm - 2*(glm::dot(wo_norm, norm))*norm);
    pdf_ret = 1.0f;
    return glm::vec3(0);
}
