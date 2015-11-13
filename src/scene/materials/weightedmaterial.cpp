#include <scene/materials/weightedmaterial.h>

WeightedMaterial::WeightedMaterial() : Material(){}
WeightedMaterial::WeightedMaterial(const glm::vec3 &color) : Material(color){}


glm::vec3 toLocalCoord(Intersection isx, glm::vec3 vec) {
    //change vectors from world->local space
    glm::mat4 worldToLocal = glm::mat4(glm::vec4(isx.tangent,0.0f),
                                       glm::vec4(isx.bitangent,0.0f),
                                       glm::vec4(isx.normal,0.0f),
                                       glm::vec4(0.0f, 0.0f,0.0f,1.0f));
    return glm::normalize(glm::vec3(glm::transpose(worldToLocal)*glm::vec4(vec,0.0f)));

}

glm::vec3 toWorldCoord(Intersection isx, glm::vec3 vec) {
    //change vectors from world->local space
    glm::mat4 localToWorld = glm::mat4(glm::vec4(isx.tangent,0.0f),
                                       glm::vec4(isx.bitangent,0.0f),
                                       glm::vec4(isx.normal,0.0f),
                                       glm::vec4(0.0f, 0.0f,0.0f,1.0f));
    return glm::normalize(glm::vec3(localToWorld*glm::vec4(vec,0.0f)));
}


glm::vec3 WeightedMaterial::EvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW,glm::vec3 &wiW_ret, float &pdf_ret, BxDFType flags) const
{
    glm::vec3 result = glm::vec3(0.0f);

    glm::vec3 local_wo = toLocalCoord(isx,woW);

    for (int i = 0; i < this->bxdfs.size(); i++) {
        //generate random number each time
        float scale=RAND_MAX+1.;
        float base=rand()/scale;
        float fine=rand()/scale;
        float rand1 = base+fine/scale;
        base = rand()/scale;
        fine = rand()/scale;
        float rand2 = base+fine/scale;

        //evaluate energy of this particular bxdf
        glm::vec3 partialEnergy = this->bxdf_weights[i]*this->bxdfs[i]->SampleAndEvaluateScatteredEnergy(local_wo, wiW_ret, rand1, rand2, pdf_ret, BSDF_ALL);
        result += partialEnergy;
    }
    wiW_ret = toWorldCoord(isx,wiW_ret);
    return result;
}

glm::vec3 WeightedMaterial::SampleAndEvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, glm::vec3 &wiW_ret, float &pdf_ret, BxDFType flags) const
{
    glm::vec3 energy = EvaluateScatteredEnergy(isx, woW, wiW_ret, pdf_ret, flags);
    pdf_ret = 1.0f;
    return energy;
}
