#include <scene/materials/bxdfs/specularreflectionBxDF.h>

glm::vec3 SpecularReflectionBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    return this->reflection_color; //perfectly specular will always reflect in direction of sampled wi!
}

float SpecularReflectionBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const{
    return 1.0f;
}

glm::vec3 SpecularReflectionBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret, BxDFType flags) const
{
    glm::vec3 norm = glm::vec3(0.0f,1.0f,0.0f);
    wi_ret = glm::normalize(wo - 2*(glm::dot(wo, norm))*norm);
    pdf_ret = 1.0f;
    return EvaluateScatteredEnergy(wo,wi_ret);
}
