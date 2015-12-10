#include "sphere.h"

#include <iostream>

#include <la.h>
#include <math.h>


static const int SPH_IDX_COUNT = 2280;  // 760 tris * 3
static const int SPH_VERT_COUNT = 382;



float UniformConePdf(float cosThetaMax)
{
    return 1.f / (2.f * PI * (1.f - cosThetaMax));
}

float Sphere::RayPDF(const Intersection &isx, const Ray &ray) {
    glm::vec3 Pcenter = transform.position();
    float radius = 0.5f*(transform.getScale().x + transform.getScale().y + transform.getScale().z)/3.0f;
    // Return uniform weight if point inside sphere
    if (glm::distance2(isx.point, Pcenter) - radius*radius < 1e-4f)
        return Geometry::RayPDF(isx, ray);

    // Compute general sphere weight
    float sinThetaMax2 = radius*radius / glm::distance2(isx.point, Pcenter);
    float cosThetaMax = glm::sqrt(glm::max(0.f, 1.f - sinThetaMax2));
    return UniformConePdf(cosThetaMax);
}


void Sphere::ComputeArea()
{
    glm::vec3 scale = transform.getScale();
    float a = scale.x*0.5f; float b = scale.y*0.5f; float c = scale.z*0.5f;
    area = 4*PI*glm::pow((glm::pow(a*b, 1.6f) + glm::pow(a*c, 1.6f) + glm::pow(b*c, 1.6f))/3.0f, 1/1.6f);
}

Intersection Sphere::SamplePoint(float a, float b, Intersection isx) {

    float z = 1.f - 2.f * a;
    float r = glm::sqrt(glm::max(0.f, 1.f - z*z));
    float phi = 2.f * PI * b;
    float x = r * glm::cos(phi);
    float y = r * glm::sin(phi);
    glm::vec3 normal3 = glm::normalize(glm::vec3(x,y,z));
    if(glm::dot(isx.normal, normal3) > 0)
    {
        normal3 = -normal3;
    }
    glm::vec4 pointL(x, y, z, 1);
    glm::vec4 normalL(normal3,0);
    glm::vec2 uv = this->GetUVCoordinates(glm::vec3(pointL));
    glm::vec3 color = Material::GetImageColor(uv, this->material->texture);
    glm::vec3 T = glm::normalize(glm::cross(glm::vec3(0,1,0), glm::vec3(normalL)));
    glm::vec3 B = glm::cross(glm::vec3(normalL), T);

    Intersection result;
    result.point = glm::vec3(transform.T() * pointL);
    result.normal = glm::normalize(glm::vec3(transform.invTransT() * normalL));
    result.texture_color = color;
    result.tangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(T, 0)));
    result.bitangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(B, 0)));
    result.object_hit = this;
    return result;

}
//    Intersection sampledPt;
//    sampledPt.normal = glm::vec3(0.0f,1.0f,0.0f);
//    glm::vec3 to_light = glm::vec3(0.0f,1.0f,0.0f);

//    while (glm::dot(to_light,sampledPt.normal) > 0) { //sample again if on backside of sphere
//        //generate another random number
//        float scale=RAND_MAX+1.;
//        float base=rand()/scale;
//        float fine=rand()/scale;
//        float rand1 = base+fine/scale;

//        float x = glm::sqrt(1-glm::pow(rand1,2))*glm::cos(2*PI*rand1);
//        float y = glm::sqrt(1-glm::pow(rand1,2))*glm::cos(2*PI*rand1);
//        float z = rand1;
//        glm::vec3 pt = glm::vec3(x,y,z); //pt on unit sphere

//        sampledPt = this->GetIntersection(Ray(glm::vec3(0.0f), glm::normalize(pt))); //intersects inside of sphere
//        //reverse direction of normal
//        sampledPt.normal = glm::normalize(sampledPt.point);
//        //offset the intersection pt
//        sampledPt.point = sampledPt.point + 0.0001f*sampledPt.normal;

//        to_light = glm::normalize(sampledPt.point - isx.point);

//}

//    //transform from local to world coords
//    sampledPt.point = glm::vec3(this->transform.T()*glm::vec4(isx.point,1.0f));
//    sampledPt.normal = glm::normalize(glm::vec3((this->transform.invTransT()*glm::vec4(isx.normal,0.0f))));
//    return sampledPt;

Intersection Sphere::GetSurfaceSample(float u1, float u2, const glm::vec3 &isx_normal)
{
    float z = 1.f - 2.f * u1;
    float r = glm::sqrt(glm::max(0.f, 1.f - z*z));
    float phi = 2.f * PI * u2;
    float x = r * glm::cos(phi);
    float y = r * glm::sin(phi);
    glm::vec3 normal3 = glm::normalize(glm::vec3(x,y,z));
    if(glm::dot(isx_normal, normal3) > 0)
    {
        normal3 = -normal3;
    }
    glm::vec4 pointL(x, y, z, 1);
    glm::vec4 normalL(normal3,0);
    glm::vec2 uv = this->GetUVCoordinates(glm::vec3(pointL));
    glm::vec3 color = Material::GetImageColor(uv, this->material->texture);
    glm::vec3 T = glm::normalize(glm::cross(glm::vec3(0,1,0), glm::vec3(normalL)));
    glm::vec3 B = glm::cross(glm::vec3(normalL), T);

    Intersection result;
    result.point = glm::vec3(transform.T() * pointL);
    result.normal = glm::normalize(glm::vec3(transform.invTransT() * normalL));
    result.texture_color = color;
    result.tangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(T, 0)));
    result.bitangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(B, 0)));
    result.object_hit = this;
    return result;
}

glm::vec3 Sphere::ComputeNormal(const glm::vec3 &P)
{}

Intersection Sphere::GetIntersection(Ray r)
{
    //Transform the ray
    Ray r_loc = r.GetTransformedCopy(transform.invT());
    Intersection result;

    float A = pow(r_loc.direction[0], 2) + pow(r_loc.direction[1], 2) + pow(r_loc.direction[2], 2);
    float B = 2*(r_loc.direction[0]*r_loc.origin[0] + r_loc.direction[1] * r_loc.origin[1] + r_loc.direction[2] * r_loc.origin[2]);
    float C = pow(r_loc.origin[0], 2) + pow(r_loc.origin[1], 2) + pow(r_loc.origin[2], 2) - 0.25f;//Radius is 0.5f
    float discriminant = B*B - 4*A*C;
    //If the discriminant is negative, then there is no real root
    if(discriminant < 0){
        return result;
    }
    float t = (-B - sqrt(discriminant))/(2*A);
    if(t < 0)
    {
        t = (-B + sqrt(discriminant))/(2*A);
    }
    if(t >= 0)
    {
        glm::vec4 P = glm::vec4(r_loc.origin + t*r_loc.direction, 1);
        result.point = glm::vec3(transform.T() * P);
        glm::vec3 normalL = glm::normalize(glm::vec3(P));
        glm::vec3 T = glm::normalize(glm::cross(glm::vec3(0,1,0), glm::vec3(normalL)));
        glm::vec3 B = glm::cross(glm::vec3(normalL), T);

        glm::vec2 uv = GetUVCoordinates(glm::vec3(P));
        result.normal = glm::normalize(glm::vec3(transform.invTransT()*glm::vec4(normalL,0.0f)));
        result.t = glm::distance(result.point, r.origin);
        result.texture_color = Material::GetImageColorInterp(uv, material->texture);
        result.object_hit = this;

        result.tangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(T, 0)));
        result.bitangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(B, 0)));
        return result;
    }
    return result;
}

glm::vec2 Sphere::GetUVCoordinates(const glm::vec3 &point)
{
    glm::vec3 p = glm::normalize(point);
    float phi = atan2f(p.z, p.x);//glm::atan(p.x/p.z);
    if(phi < 0)
    {
        phi += TWO_PI;
    }
    float theta = glm::acos(p.y);
    return glm::vec2(1 - phi/TWO_PI, 1 - theta / PI);
}

// These are functions that are only defined in this cpp file. They're used for organizational purposes
// when filling the arrays used to hold the vertex and index data.
void createSphereVertexPositions(glm::vec3 (&sph_vert_pos)[SPH_VERT_COUNT])
{
    // Create rings of vertices for the non-pole vertices
    // These will fill indices 1 - 380. Indices 0 and 381 will be filled by the two pole vertices.
    glm::vec4 v;
    // i is the Z axis rotation
    for (int i = 1; i < 19; i++) {
        // j is the Y axis rotation
        for (int j = 0; j < 20; j++) {
            v = glm::rotate(glm::mat4(1.0f), j / 20.f * TWO_PI, glm::vec3(0, 1, 0))
                    * glm::rotate(glm::mat4(1.0f), -i / 18.0f * PI, glm::vec3(0, 0, 1))
                    * glm::vec4(0, 0.5f, 0, 1);
            sph_vert_pos[(i - 1) * 20 + j + 1] = glm::vec3(v);
        }
    }
    // Add the pole vertices
    sph_vert_pos[0] = glm::vec3(0, 0.5f, 0);
    sph_vert_pos[381] = glm::vec3(0, -0.5f, 0);  // 361 - 380 are the vertices for the bottom cap
}


void createSphereVertexNormals(glm::vec3 (&sph_vert_nor)[SPH_VERT_COUNT])
{
    // Unlike a cylinder, a sphere only needs to be one normal per vertex
    // because a sphere does not have sharp edges.
    glm::vec4 v;
    // i is the Z axis rotation
    for (int i = 1; i < 19; i++) {
        // j is the Y axis rotation
        for (int j = 0; j < 20; j++) {
            v = glm::rotate(glm::mat4(1.0f), j / 20.0f * TWO_PI, glm::vec3(0, 1, 0))
                    * glm::rotate(glm::mat4(1.0f), -i / 18.0f * PI, glm::vec3(0, 0, 1))
                    * glm::vec4(0, 1.0f, 0, 0);
            sph_vert_nor[(i - 1) * 20 + j + 1] = glm::vec3(v);
        }
    }
    // Add the pole normals
    sph_vert_nor[0] = glm::vec3(0, 1.0f, 0);
    sph_vert_nor[381] = glm::vec3(0, -1.0f, 0);
}


void createSphereIndices(GLuint (&sph_idx)[SPH_IDX_COUNT])
{
    int index = 0;
    // Build indices for the top cap (20 tris, indices 0 - 60, up to vertex 20)
    for (int i = 0; i < 19; i++) {
        sph_idx[index] = 0;
        sph_idx[index + 1] = i + 1;
        sph_idx[index + 2] = i + 2;
        index += 3;
    }
    // Must create the last triangle separately because its indices loop
    sph_idx[57] = 0;
    sph_idx[58] = 20;
    sph_idx[59] = 1;
    index += 3;

    // Build indices for the body vertices
    // i is the Z axis rotation
    for (int i = 1; i < 19; i++) {
        // j is the Y axis rotation
        for (int j = 0; j < 20; j++) {
            sph_idx[index] = (i - 1) * 20 + j + 1;
            sph_idx[index + 1] = (i - 1) * 20 +  j + 2;
            sph_idx[index + 2] = (i - 1) * 20 +  j + 22;
            sph_idx[index + 3] = (i - 1) * 20 +  j + 1;
            sph_idx[index + 4] = (i - 1) * 20 +  j + 22;
            sph_idx[index + 5] = (i - 1) * 20 +  j + 21;
            index += 6;
        }
    }

    // Build indices for the bottom cap (20 tris, indices 2220 - 2279)
    for (int i = 0; i < 19; i++) {
        sph_idx[index] = 381;
        sph_idx[index + 1] = i + 361;
        sph_idx[index + 2] = i + 362;
        index += 3;
    }
    // Must create the last triangle separately because its indices loop
    sph_idx[2277] = 381;
    sph_idx[2278] = 380;
    sph_idx[2279] = 361;
    index += 3;
}

void Sphere::create()
{
    GLuint sph_idx[SPH_IDX_COUNT];
    glm::vec3 sph_vert_pos[SPH_VERT_COUNT];
    glm::vec3 sph_vert_nor[SPH_VERT_COUNT];
    glm::vec3 sph_vert_col[SPH_VERT_COUNT];

    createSphereVertexPositions(sph_vert_pos);
    createSphereVertexNormals(sph_vert_nor);
    createSphereIndices(sph_idx);
    for (int i = 0; i < SPH_VERT_COUNT; i++) {
        sph_vert_col[i] = material->base_color;
    }

    count = SPH_IDX_COUNT;

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(sph_idx, SPH_IDX_COUNT * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(sph_vert_pos, SPH_VERT_COUNT * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(sph_vert_col, SPH_VERT_COUNT * sizeof(glm::vec3));

    bufNor.create();
    bufNor.bind();
    bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufNor.allocate(sph_vert_nor, SPH_VERT_COUNT * sizeof(glm::vec3));

    this->boundingBox = new BoundingBox();
    setBoundingBox();
}

/**
 * @brief Cube::setBoundingBox
 * @brief called from after this.transform is defined
 * @brief and TRANSFORMED (max x, max y, max z)
 *
 * @param cub_vert_pos (vertex positions of cubes)
 */
void Sphere::setBoundingBox() {
    this->boundingBox->setTransformedBox(this->transform.T());
    boundingBox->create();
}


bool Sphere::isMesh() {
    return false;
}
