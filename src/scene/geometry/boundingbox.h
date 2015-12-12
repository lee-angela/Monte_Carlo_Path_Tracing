#ifndef BOUNDINGBOX
#define BOUNDINGBOX

#include <openGL/drawable.h>
#include <raytracing/intersection.h>

class Intersection;

class BoundingBox : public Drawable
{

public:


    GLuint box_idx[24];
    glm::vec3 box_vert_pos[8];
    glm::vec3 box_vert_col[8];


    glm::mat4 transformation; //transformation matrices of this world axis-aligned bounding box
    glm::mat4 invTransformation;
    glm::mat4 invTransposeT;

    BoundingBox();
    BoundingBox(glm::vec3 min, glm::vec3 max);
    ~BoundingBox();

    //coordinates stored in world space
    glm::vec3 minBound; //min bound of geom (or BVH node)
    glm::vec3 maxBound; //max bound of "
    glm::vec3 centerPoint; //avg of bounds

    void setTransformedBox(glm::mat4 transformation); //sets bbox bounds for geometries, using transformation matrix (object-aligned)
    void setWorldBox(); //sets bounds for world axis aligned box
    void combineBoxes(BoundingBox* a); //returns new combination box using a and b's min/max bounds

    void createBoxVertexPositions();
    void createBoxIndices();

    glm::vec2 GetUVCoordinates(const glm::vec3 &point);
    void setBoundingBox(glm::mat4 transf);
    void setTransformationMats(glm::mat4 scale, glm::mat4 transl);

    Intersection getIntersection(Ray r);
    float getSurfArea();

    void create();
    virtual GLenum drawMode();
    bool isMesh();

};

#endif // BOUNDINGBOX

