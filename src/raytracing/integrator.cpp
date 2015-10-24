#include <raytracing/integrator.h>
#include <math.h>


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

//Basic ray trace
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{
    if (depth >= max_depth) { //do not continue if max depth is reached
        return glm::vec3(0.0f);
    }

    glm::vec3 light_transport;
    //light_transport = emitted_radiance + brdf(p, ray_in, ray_out) * incoming_light * dot(incoming_light,p.N)
    glm::vec3 rad_emit;
    glm::vec3 brdf;
    Ray incoming_ray;
    glm::vec3 incoming_light;

    Intersection pt = this->intersection_engine->GetIntersection(r); //find intersection from ray


    //pick random light source to trace ray to
    int random = floor(rand()%this->scene->lights.count());
    //pick random pt on this light's geometry to trace ray to
    this->scene->lights[random]->SamplePoint(rand(), rand());


    if (pt.object_hit->material->is_light_source) {
        rad_emit = pt.object_hit->material->EvaluateScatteredEnergy(pt,glm::vec3(0.0f), glm::vec3(0.0f),BSDF_ALL);
    }

    light_transport = glm::vec3(0.0f);




    return glm::vec3(0,0,0);
}

void Integrator::SetDepth(unsigned int depth)
{
    max_depth = depth;
}
