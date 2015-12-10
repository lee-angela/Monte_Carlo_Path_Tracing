#include <scene/geometry/disc.h>

void Disc::ComputeArea()
{
    //r=0.5 of unit disk

    //A = PI*a*b (if a=b=r, then PI*r*r)

    //find transformed r:
    glm::mat4 T = this->transform.T();
    float a = glm::length(T*glm::vec4(0.5f,0.0f,0.0f,1.0f) - T*glm::vec4(0.0f,0.0f,0.0f,1.0f));
    float b = glm::length(T*glm::vec4(0.0f,0.5f,0.0f,1.0f) - T*glm::vec4(0.0f,0.0f,0.0f,1.0f));

    this->area = PI*a*b;
}


/**
 * @brief SamplePoint
 * @param float a (rand float between 0 and 1)
 * @param float b (rand float between 0 and 1)
 * @return Intersection in world coordinates of the sampled point
 */
Intersection Disc::SamplePoint(float a, float b, Intersection isx) {
    glm::vec4 pt1 = glm::vec4(0.5f,0.0f,0.0f,1.0f); //any point on circumfrence
    float angle = a*2*PI; //fraction of the entire circumfrence of disc
    glm::vec3 axis(0,0,1); //rotate around z axis

    glm::vec3 pt2 = glm::vec3(glm::rotate(glm::mat4(1.0f), angle, axis) * pt1);
    glm::vec4 finalPt;
    finalPt = glm::vec4(b*pt2,1.0f); //between pt2 and origin
    finalPt = this->transform.T()*finalPt; //transform this sampled pt

    Intersection pt;
    pt.point = glm::vec3(finalPt);
    pt.normal = glm::vec3(glm::normalize(this->transform.invTransT()*glm::vec4(0.0f,0.0f,1.0f,0.0f)));
    pt.object_hit = this;
    return pt;
}



Intersection Disc::GetIntersection(Ray r)
{
    //Transform the ray
    Ray r_loc = r.GetTransformedCopy(transform.invT());
    Intersection result;

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    glm::vec4 P = glm::vec4(t * r_loc.direction + r_loc.origin, 1);
    //Check that P is within the bounds of the disc (not bothering to take the sqrt of the dist b/c we know the radius)
    float dist2 = (P.x * P.x + P.y * P.y);
    if(t > 0 && dist2 <= 0.25f)
    {
        glm::vec3 normalL = glm::vec3(0.0f,0.0f,1.0f);
        result.point = glm::vec3(transform.T() * P);
        result.normal = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(ComputeNormal(glm::vec3(P)), 0)));
        result.object_hit = this;
        result.t = glm::distance(result.point, r.origin);
        result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(glm::vec3(P)), material->texture);
        //TODO: Store the tangent and bitangent
        glm::vec3 T = glm::normalize(glm::cross(glm::vec3(0,1,0), glm::vec3(normalL)));
        glm::vec3 B = glm::cross(glm::vec3(normalL), T);

        result.normal = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(normalL, 0.0f)));
        result.tangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(T, 0)));
        result.bitangent = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(B, 0)));
        return result;
    }
    return result;
}

glm::vec2 Disc::GetUVCoordinates(const glm::vec3 &point)
{
    return glm::vec2(point.x + 0.5f, point.y + 0.5f);
}

glm::vec3 Disc::ComputeNormal(const glm::vec3 &P)
{
    return glm::vec3(0,0,1);
}

void Disc::create()
{
    GLuint idx[54];
    //18 tris, 54 indices
    glm::vec3 vert_pos[20];
    glm::vec3 vert_nor[20];
    glm::vec3 vert_col[20];

    //Fill the positions, normals, and colors
    glm::vec4 pt(0.5f, 0, 0, 1);
    float angle = 18.0f * DEG2RAD;
    glm::vec3 axis(0,0,1);
    for(int i = 0; i < 20; i++)
    {
        //Position
        glm::vec3 new_pt = glm::vec3(glm::rotate(glm::mat4(1.0f), angle * i, axis) * pt);
        vert_pos[i] = new_pt;
        //Normal
        vert_nor[i] = glm::vec3(0,0,1);
        //Color
        vert_col[i] = material->base_color;
    }

    //Fill the indices.
    int index = 0;
    for(int i = 0; i < 18; i++)
    {
        idx[index++] = 0;
        idx[index++] = i + 1;
        idx[index++] = i + 2;
    }

    count = 54;

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(idx, 54 * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(vert_pos, 20 * sizeof(glm::vec3));

    bufNor.create();
    bufNor.bind();
    bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufNor.allocate(vert_nor, 20 * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(vert_col, 20 * sizeof(glm::vec3));
}

/**
 * @brief Cube::setBoundingBox
 * @brief called from after this.transform is defined
 * @brief and TRANSFORMED (max x, max y, max z)
 *
 * @param cub_vert_pos (vertex positions of cubes)
 */
void Disc::setBoundingBox() {
    this->boundingBox->setTransformedBox(this->transform.T());
    boundingBox->create();
}


bool Disc::isMesh() {
    return false;
}
