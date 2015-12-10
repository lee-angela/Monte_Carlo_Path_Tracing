#include <scene/geometry/square.h>

void SquarePlane::ComputeArea()
{
    glm::vec4 side1 = this->transform.T()*glm::vec4(-0.5f,-0.5f,0.0f,1.0f) -
                        this->transform.T()*glm::vec4(0.5f,-0.5f,0.0f,1.0f); //left edge
    glm::vec4 side2 = this->transform.T()*glm::vec4(0.5f,-0.5f,0.0f,1.0f) -
                        this->transform.T()*glm::vec4(0.5f,0.5f,0.0f,1.0f); //front edge
    area = glm::length(glm::cross(glm::vec3(side1), glm::vec3(side2)));
}

/**
 * @brief SamplePoint
 * @param float a (rand float between 0 and 1)
 * @param float b (rand float between 0 and 1)
 * @return Intersection in world coordinates of the sampled point
 */
Intersection SquarePlane::SamplePoint(float a, float b, Intersection isx) {
    Intersection pt;

    //start along one edge
    glm::vec3 pt1 = glm::vec3(-0.5+b,-0.5+a,0);
    pt.point = glm::vec3(this->transform.T()*glm::vec4(pt1, 1.0f));
    pt.normal = glm::normalize(glm::vec3(this->transform.invTransT()*glm::vec4(0.0f,0.0f,1.0f,0.0f)));
    pt.object_hit = this;

    return pt;
}

Intersection SquarePlane::GetIntersection(Ray r)
{
    //Transform the ray
    Ray r_loc = r.GetTransformedCopy(transform.invT());
    Intersection result;

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    glm::vec4 P = glm::vec4(t * r_loc.direction + r_loc.origin, 1);
    //Check that P is within the bounds of the square
    if(t > 0 && P.x >= -0.5f && P.x <= 0.5f && P.y >= -0.5f && P.y <= 0.5f)
    {
        glm::vec3 normalL = glm::vec3(0.0f,0.0f,1.0f);
        result.point = glm::vec3(transform.T() * P);
        result.normal = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(normalL, 0.0f)));
        result.object_hit = this;
        result.t = glm::distance(result.point, r.origin);
        result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(glm::vec3(P)), material->texture);
        //TODO: Store the tangent and bitangent
        glm::vec3 T = glm::normalize(glm::cross(glm::vec3(0.0f,1.0f,0.0f), normalL));
        glm::vec3 B = glm::normalize(glm::cross(glm::vec3(normalL), T));
        result.tangent = glm::normalize(glm::vec3(transform.invTransT()*glm::vec4(T,0.0f)));
        result.bitangent = glm::normalize(glm::vec3(transform.invTransT()*glm::vec4(B, 0.0f)));

        return result;
    }
    return result;
}


glm::vec2 SquarePlane::GetUVCoordinates(const glm::vec3 &point)
{
    return glm::vec2(point.x + 0.5f, point.y + 0.5f);
}

glm::vec3 SquarePlane::ComputeNormal(const glm::vec3 &P)
{
        return glm::vec3(0,0,1);
}

void SquarePlane::create()
{
    GLuint cub_idx[6];
    glm::vec3 cub_vert_pos[4];
    glm::vec3 cub_vert_nor[4];
    glm::vec3 cub_vert_col[4];

    cub_vert_pos[0] = glm::vec3(-0.5f, 0.5f, 0);  cub_vert_nor[0] = glm::vec3(0, 0, 1); cub_vert_col[0] = material->base_color;
    cub_vert_pos[1] = glm::vec3(-0.5f, -0.5f, 0); cub_vert_nor[1] = glm::vec3(0, 0, 1); cub_vert_col[1] = material->base_color;
    cub_vert_pos[2] = glm::vec3(0.5f, -0.5f, 0);  cub_vert_nor[2] = glm::vec3(0, 0, 1); cub_vert_col[2] = material->base_color;
    cub_vert_pos[3] = glm::vec3(0.5f, 0.5f, 0);   cub_vert_nor[3] = glm::vec3(0, 0, 1); cub_vert_col[3] = material->base_color;

    cub_idx[0] = 0; cub_idx[1] = 1; cub_idx[2] = 2;
    cub_idx[3] = 0; cub_idx[4] = 2; cub_idx[5] = 3;

    count = 6;

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(cub_idx, 6 * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(cub_vert_pos, 4 * sizeof(glm::vec3));

    bufNor.create();
    bufNor.bind();
    bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufNor.allocate(cub_vert_nor, 4 * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(cub_vert_col, 4 * sizeof(glm::vec3));

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
void SquarePlane::setBoundingBox() {
    this->boundingBox->setTransformedBox(this->transform.T());
    boundingBox->create();
}


bool SquarePlane::isMesh() {
    return false;
}
