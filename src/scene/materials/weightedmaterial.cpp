#include <scene/materials/weightedmaterial.h>

WeightedMaterial::WeightedMaterial() : Material(){}
WeightedMaterial::WeightedMaterial(const glm::vec3 &color) : Material(color){}

glm::vec3 WeightedMaterial::EvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, const glm::vec3 &wiW, BxDFType flags) const
{
    glm::vec3 result = glm::vec3(0.0f);
    float pdf_ret;

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
        glm::vec3 partialEnergy = this->bxdf_weights[i]*this->bxdfs[i]->SampleAndEvaluateScatteredEnergy(woW, wiW, rand1, rand2, pdf_ret);
        result += partialEnergy;
    }
    return result;
}

glm::vec3 WeightedMaterial::SampleAndEvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, glm::vec3 &wiW_ret, float &pdf_ret, BxDFType flags) const
{
    glm::vec3 energy = EvaluateScatteredEnergy(isx, woW, wiW_ret, flags);
    pdf_ret = 1.0f;
    return energy;
}
