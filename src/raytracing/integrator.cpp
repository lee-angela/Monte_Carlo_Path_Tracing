#include <raytracing/integrator.h>
#include <math.h>
#include <chrono>


Integrator::Integrator():
    max_depth(5)
{
    scene = NULL;
    intersection_engine = NULL;
}

glm::vec3 ComponentMult(const glm::vec3 &a, const glm::vec3 &b)
{
    return glm::vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}


float generate_random() {
    float scale=RAND_MAX+1.;
    float base=rand()/scale;
    float fine=rand()/scale;
    return (base+fine)/scale;
}


glm::vec3 toLocalCoords(Intersection isx, glm::vec3 vec) {
    //change vectors from world->local space
    glm::mat4 worldToLocal = glm::mat4(glm::vec4(isx.tangent,0.0f),
                                       glm::vec4(isx.bitangent,0.0f),
                                       glm::vec4(isx.normal,0.0f),
                                       glm::vec4(0.0f, 0.0f,0.0f,1.0f));
    return glm::normalize(glm::vec3(glm::transpose(worldToLocal)*glm::vec4(vec,0.0f)));

}

glm::vec3 toWorldCoords(Intersection isx, glm::vec3 vec) {
    //change vectors from world->local space
    glm::mat4 localToWorld = glm::mat4(glm::vec4(isx.tangent,0.0f),
                                       glm::vec4(isx.bitangent,0.0f),
                                       glm::vec4(isx.normal,0.0f),
                                       glm::vec4(0.0f, 0.0f,0.0f,1.0f));
    return glm::normalize(glm::vec3(localToWorld*glm::vec4(vec,0.0f)));
}


glm::vec3 Integrator::EstimateDirectLighting(const Intersection &isx, unsigned int &samples_taken) {

    if (isx.object_hit == NULL) {
        return glm::vec3(0.0f);
    }

  //CALCULATE LIGHT MIS (first half of total direct lighting LTE)
    glm::vec3 LD_light_sampling;
    glm::vec3 wo = glm::normalize(this->scene->camera.eye - isx.point);


    //pick random light source to trace ray to
    int random = floor(rand()%this->scene->lights.count());
    Intersection lightSample = this->scene->lights[random]->SamplePoint(generate_random(),generate_random(),isx);

    //incoming ray of light (towards light, originating from isx.point)
    glm::vec3 wi = glm::normalize(lightSample.point - isx.point);

     int bxdf_idx;
    //check if there are obstructions
    Ray light_feeler = Ray(isx.point, glm::normalize(wi));
    //if the object hit is the light source, return light's energy calc
    if (isx.object_hit->material->is_light_source) {
        glm::vec3 rad_emit = isx.object_hit->material->EvaluateScatteredEnergy(isx,-wi, -wi,bxdf_idx,BSDF_ALL);
        return rad_emit;
    }

    glm::vec3 brdf = isx.object_hit->material->EvaluateScatteredEnergy(isx,toLocalCoords(isx,wo),toLocalCoords(isx,wi),bxdf_idx,BSDF_ALL);
    float bxdf_pdf = isx.object_hit->material->bxdfs[bxdf_idx]->PDF(toLocalCoords(isx,wo),toLocalCoords(isx,wi));

    //flip the direction of light_feeler ray that is used to calculate the pdf of light
    Ray pdfLightRay = Ray(isx.point, glm::normalize(lightSample.point-isx.point));

    //isx is on light
    float light_pdf = lightSample.object_hit->RayPDF(lightSample, pdfLightRay);
    //find power heuristic
    float light_pwr_heuristic = glm::pow(light_pdf, 2.0f)/ (glm::pow(light_pdf, 2.0f) + glm::pow(bxdf_pdf,2.0f));

    //test for shadows
    Intersection obstruction;
    for (int i = 0; i < this->scene->objects.size(); i++) {
        Intersection temp = this->scene->objects[i]->GetIntersection(light_feeler);
        if ((temp.t < obstruction.t && temp.t > 0) || (temp.t > 0 && obstruction.t < 0)) {
            obstruction = temp;
        }
    }
    //if this light feeler hits a different object, return black color
    glm::vec3 Li;
    if (obstruction.object_hit != lightSample.object_hit) {
        Li = glm::vec3(0.0f,0.0f,0.0f);
    } else {
        Li = lightSample.object_hit->material->EvaluateScatteredEnergy(lightSample,-wi,-wi,bxdf_idx,BSDF_ALL);
    }

    if (light_pdf < 0.0001) {
        LD_light_sampling = glm::vec3(0.0f);
    } else if (light_pwr_heuristic > 0.0001f) {
        float temp = glm::abs(glm::dot(toWorldCoords(isx,wi), isx.normal));
        LD_light_sampling = brdf * Li * glm::abs(glm::dot(wi, isx.normal)) * light_pwr_heuristic /light_pdf;
    } else {
        LD_light_sampling = glm::vec3(0.0f);
    }



  //CALCULATE BRDF MIS (second half of direct lighting LTE)
    glm::vec3 LD_brdf_sampling;
    glm::vec3 bxdf = isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx, toLocalCoords(isx,wo), wi, bxdf_pdf ,BSDF_ALL);
    wo = toWorldCoords(isx,wo);

    //check if the wi found (towards light) intersects a light source
    Intersection brdf_sampled_light;
    for (int i = 0; i < this->scene->objects.size(); i++) {
        Intersection temp = this->scene->objects[i]->GetIntersection(light_feeler);
        if ((temp.t < obstruction.t && temp.t > 0) || (temp.t > 0 && obstruction.t < 0)) {
            brdf_sampled_light = temp;
        }
    }
    // if sampled wi does not intersect light source (i.e. there is obstruction or never reaches light)
    // this is black
    if (brdf_sampled_light.object_hit == NULL || !brdf_sampled_light.object_hit->material->is_light_source) {
        LD_brdf_sampling = glm::vec3(0.0f);
    } else {
        wi = glm::normalize(brdf_sampled_light.point - isx.point);
        Li = brdf_sampled_light.object_hit->material->EvaluateScatteredEnergy(
                    brdf_sampled_light,toLocalCoords(brdf_sampled_light,wo),toLocalCoords(brdf_sampled_light,wi),bxdf_idx,BSDF_ALL);

        float bxdf_pwr_heuristic = glm::pow(bxdf_pdf, 2)/ (glm::pow(light_pdf, 2) + glm::pow(bxdf_pdf,2));

        LD_brdf_sampling = bxdf*Li*glm::abs(glm::dot(isx.normal,wi))*bxdf_pwr_heuristic/bxdf_pdf;
    }


    return LD_light_sampling + LD_brdf_sampling;
}


glm::vec3 Integrator::EstimateIndirectLighting(const Intersection &isx, unsigned int &samples_taken) {
    // specular reflective - will not be seen by direct lighting because there is 0 probability that the sampled ray
    // (i.e. the reflection ray) will see the randomly sampled point on the light

    float throughput = 1.0f;
    int depth = 0;
    glm::vec3 ID_lighting;

    if (isx.object_hit->material->is_light_source) {
        return glm::vec3(0.0f); //light emitted is calculated in EstimateDirectLighting()!
    }
    //if wi sampled sees a light, ignore -> 0 (at depth 0)
    if (depth == 0) {
        //sample wi and see if it sees a light
        bool seesLight = false;
        glm::vec3 woW = glm::normalize(this->scene->camera.eye - isx.point);
        glm::vec3 sampled_wi;
        float pdf_ret;

        //transform the wo and wi to local space before sending to evaluateEnergy functions
        isx.object_hit->material->SampleAndEvaluateScatteredEnergy(isx,woW,sampled_wi, pdf_ret,BSDF_ALL);

        for (int i = 0; i < this->scene->lights.size(); i++) {
            Intersection light_isx = this->scene->lights[i]->GetIntersection(Ray(isx.point, sampled_wi));
            if (light_isx.object_hit != NULL) {
                seesLight = true;
            }
        }
        if (seesLight) {
            return glm::vec3(0.0f);
        }
    }

    //light estimation loop
    //during iteration, keep track of all comps of LTE except for the Li (light) (energy is multiplicative)
    //this calculates all the things that are touched/light carried throughout scene
    //at end of iterations, if hits light, Li = light.ESE (just multiply into the rest of LTE eq that was kept track of)
    // if hits nothing, Li = 0, ends in darkness
    while (true) {
        depth++;
        if (depth >= max_depth) {
            break;
        }

        //LTE - calculate the equation
        //calculate the brdf()*|dot| component of the LTE (light transport equation)
//        glm::vec3 rgb_comp = brdf()*glm::abs(glm::dot());
//        glm::vec3 Li = lightSample.object_hit->material->EvaluateScatteredEnergy(lightSample,-wi,-wi,BSDF_ALL); //light radiance

//        float maxRGB = glm::max(rgb_comp[2],glm::max(rgb_comp[0], rbg_comp[1])); //find max comp of this vector



        //Russian Roulette
        if (depth > 2) {
            float comparator = generate_random(); //get uniform rand number
            if (comparator > throughput) {
                break; //if rand number > throughput, stop ray.
            }
        }

        //if this iteration NOT filtered with Russian Roulette, then continue to find the indirect lighting

        //redefined sampled_wi to use in next iteration
//        sampled_wi = ;


    }
    return ID_lighting;

}



//Basic ray trace
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{
    if (depth >= max_depth) { //do not continue if max depth is reached
        return glm::vec3(0.0f);
    }
    Intersection isx = this->intersection_engine->GetIntersection(r); //find intersection from ray
    isx.point = isx.point + 0.0001f*isx.normal;

    unsigned int N = 10;
    glm::vec3 direct_light = EstimateDirectLighting(isx,N);
//    glm::vec3 indirect_light = EstimateIndirectLighting(isx,N);

    return direct_light/*+indirect_light*/;
}



void Integrator::SetDepth(unsigned int depth)
{
    max_depth = depth;
}
