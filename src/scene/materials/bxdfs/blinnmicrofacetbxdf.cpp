#include <scene/materials/bxdfs/blinnmicrofacetbxdf.h>
#include <math.h>

glm::vec3 BlinnMicrofacetBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //microfacets act as specular reflective
    //normal (or H direction, N of the microfacet) is in local space as y axis
    return this->reflection_color;
}

float BlinnMicrofacetBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const {
    //wo and wi must be in local space

    float scale=RAND_MAX+1.;
    float base=rand()/scale;
    float fine=rand()/scale;
    float rand1= (base+fine)/scale;
    base=rand()/scale;
    fine=rand()/scale;
    float rand2 = (base+fine)/scale;

    //sample half-angle vector wh
    float costheta = glm::pow(rand1, 1/(exponent+1));
    float sintheta = glm::sqrt(glm::max(0.0f, 1.0f-costheta*costheta));
    float phi = rand2*2.0f*PI;
    glm::vec3 wh = glm::vec3(sintheta*glm::cos(phi), sintheta*glm::sin(phi), costheta);

    if (!(wo[2]*wh[2] > 0.0f)) {
        wh = -wh;
    }
    //compute pdf for wi from blinn distribution
   float pdf_ret = ((exponent+1.0f)*glm::pow(costheta, exponent)) / (2.0f*PI*4.0f*glm::dot(wo,wh));
    if (glm::dot(wo, wh) <= 0.0f) {
        pdf_ret = 0.0f;
    }
    return pdf_ret;
}

glm::vec3 BlinnMicrofacetBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret,float rand1, float rand2, float &pdf_ret, BxDFType flags) const
{ // use Torrance-Sparrow model

    //sample half-angle vector wh
    float costheta = glm::pow(rand1, 1/(exponent+1));
    float sintheta = glm::sqrt(glm::max(0.0f, 1.0f-costheta*costheta));
    float phi = rand2*2.0f*PI;
    glm::vec3 wh = glm::vec3(sintheta*glm::cos(phi), sintheta*glm::sin(phi), costheta);

    if (!(wo[2]*wh[2] > 0.0f)) {
        wh = -wh;
    }

    //compute incident direction by reflecting about wh
    wi_ret = -wo + 2.0f*glm::dot(wo, wh)*wh; //currently in local space

    //compute pdf for wi from blinn distribution
    pdf_ret = ((exponent+1.0f)*glm::pow(costheta, exponent)) / (2.0f*PI*4.0f*glm::dot(wo,wh));
    if (glm::dot(wo, wh) <= 0.0f) {
        pdf_ret = 0.0f;
    }

    return EvaluateScatteredEnergy(wo, wi_ret);
}
