#include <scene/materials/bxdfs/blinnmicrofacetbxdf.h>
#include <math.h>

glm::vec3 BlinnMicrofacetBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //microfacets act as specular reflective
    //normal (or H direction, N of the microfacet) is in local space as y axis
    return this->reflection_color;
}


glm::vec3 BxDF::SampleAndEvaluateScatteredEnergy(Intersection isx, const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{ // use Torrance-Sparrow model

    //change vectors from world->local space
    glm::mat4 worldToLocal = glm::mat4(glm::vec4(isx.tangent,0.0f),
                                    glm::vec4(isx.bitangent,0.0f),
                                    glm::vec4(isx.normal,0.0f),
                                    glm::vec4(0.0f, 0.0f,0.0f,1.0f));
    glm::vec3 local_wo = glm::vec3(worldToLocal*glm::vec4(vec,0.0f));


    glm::vec3 norm = glm::vec3(0.0f,1.0f,0.0f);
    //sample direction H that microfacet is facing
    glm::vec3 H = glm::vec3(0.0f);

    //change wo and wi to make H the local y axis
    glm::mat4 transf = glm::mat4();
    wi_ret = glm::vec3(transf*glm::vec4(0.0f));
    glm::vec3 energy = EvaluateScatteredEnergy(wo,wi_ret);

    //calculate distribution that microfacet will be in H direction
    float dotted = glm::dot(H, norm);
    float D = (exponent + 2) / (2*PI*glm::pow(dotted, exponent));

    wi_ret = glm::vec3(0.0f); //currently in local space
    //transform the wi_ret back to world space
    glm::mat4 localToWorld = glm::transpose(worldToLocal);
    wi_ret = localToWorld*wi_ret;

    pdf_ret = ;
    return glm::vec3(0);
}
