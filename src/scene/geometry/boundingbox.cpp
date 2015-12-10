#include <scene/geometry/BoundingBox.h>
#include <scene/geometry/geometry.h>
#include <raytracing/intersection.h>


#include <iostream>
#include <sstream>
using namespace std;


BoundingBox::BoundingBox() :
    minBound(glm::vec3(0.0f)),
    maxBound(glm::vec3(0.0f)),
    centerPoint(glm::vec3(0.0f))
{
    this->create();
}


BoundingBox::~BoundingBox() {
}

void BoundingBox::setTransformationMats(glm::mat4 scale, glm::mat4 transl) {
    this->transformation = transl*scale;
    this->invTransformation = glm::inverse(transformation);
    this->invTransposeT = glm::inverse(glm::transpose(this->transformation));
}



/**
 * @brief BoundingBox::setTransformedBox
 * @brief transforms the box by the transformation matrix input
 * @brief sets the min/max bounds of world aligned bounding box to that of the
 * @param transformation
 */
void BoundingBox::setTransformedBox(glm::mat4 transformation) {
    float Xminmax[2];
    float Yminmax[2];
    float Zminmax[2];

    Xminmax[0] = (transformation*glm::vec4(box_vert_pos[0],1.0f))[0];
    Xminmax[1] = (transformation*glm::vec4(box_vert_pos[0],1.0f))[0];

    Yminmax[0] = (transformation*glm::vec4(box_vert_pos[0],1.0f))[1];
    Yminmax[1] = (transformation*glm::vec4(box_vert_pos[0],1.0f))[1];

    Zminmax[0] = (transformation*glm::vec4(box_vert_pos[0],1.0f))[2];
    Zminmax[1] = (transformation*glm::vec4(box_vert_pos[0],1.0f))[2];

    for (int i = 0; i < 8; i++) {
        glm::vec3 newVrt = glm::vec3(transformation*glm::vec4(box_vert_pos[i],1.0f));
        if (newVrt[0] < Xminmax[0]) {
            Xminmax[0] = newVrt[0];
        }
        if (newVrt[0] > Xminmax[1]) {
            Xminmax[1] = newVrt[0];
        }
        if (newVrt[1] < Yminmax[0]) {
            Yminmax[0] = newVrt[1];
        }
        if (newVrt[1] > Yminmax[1]) {
            Yminmax[1] = newVrt[1];
        }
        if (newVrt[2] < Zminmax[0]) {
            Zminmax[0] = newVrt[2];
        }
        if (newVrt[2] > Zminmax[1]) {
            Zminmax[1] = newVrt[2];
        }
    } //end for loop through all box vertices

    this->minBound = glm::vec3(Xminmax[0], Yminmax[0], Zminmax[0]);
    this->maxBound = glm::vec3(Xminmax[1], Yminmax[1], Zminmax[1]);
    this->centerPoint = glm::vec3((this->minBound[0]+this->maxBound[0])/2, (this->minBound[1]+this->maxBound[1])/2,
            (this->minBound[2]+this->maxBound[2])/2);

    //find new transformation
    float Xscale = glm::abs(Xminmax[1] - Xminmax[0]); // width of box in any dimension
    float Yscale = glm::abs(Yminmax[1] - Yminmax[0]);
    float Zscale = glm::abs(Zminmax[1] - Zminmax[0]);
    setTransformationMats(glm::scale(glm::mat4(1.0f),glm::vec3(Xscale, Yscale, Zscale)), glm::translate(glm::mat4(1.0f), this->centerPoint));

}





/**
 * @brief BoundingBox::combineBoxes
 * @brief makes this bounding box a larger bounding box that is the combination of a and b bounding boxes.
 * @param a
 * @param b
 */
void BoundingBox::combineBoxes(BoundingBox* a) {
    float Xmin = glm::min(a->minBound[0], this->minBound[0]);
    float Xmax = glm::max(a->maxBound[0], this->maxBound[0]);

    float Ymin = glm::min(a->minBound[1], this->minBound[1]);
    float Ymax = glm::max(a->maxBound[1], this->maxBound[1]);

    float Zmin = glm::min(a->minBound[2], this->minBound[2]);
    float Zmax = glm::max(a->maxBound[2], this->maxBound[2]);


    //recalculate min/max/centerpoint of this combined bounding box
    this->minBound = glm::vec3(Xmin, Ymin, Zmin);
    this->maxBound = glm::vec3(Xmax, Ymax, Zmax);
    this->centerPoint = glm::vec3((this->minBound[0]+this->maxBound[0])/2, (this->minBound[1]+this->maxBound[1])/2,
            (this->minBound[2]+this->maxBound[2])/2);

    //find new transformation
    float Xscale = glm::abs(this->maxBound[0] - this->minBound[0]); //width of box in any dimension
    float Yscale = glm::abs(this->maxBound[1] - this->minBound[1]);
    float Zscale = glm::abs(this->maxBound[2] - this->minBound[2]);

    //set transformation matrix of this combined box
    setTransformationMats(glm::scale(glm::mat4(1.0f),glm::vec3(Xscale, Yscale, Zscale)), glm::translate(glm::mat4(1.0f), this->centerPoint));

}


Intersection BoundingBox::getIntersection(Ray r) {

    // bounding boxes will be world aligned, with min/max bounds defined

    Ray r1 = r.GetTransformedCopy(this->invTransformation);
    //returns the closest point of intersection from the ray's origin

    float t_near = -INFINITY;
    float t_far = INFINITY;
    //default Intersection pt
    Intersection pt;
    glm::vec4 normNear;
    glm::vec4 normFar;

    float t0;
    float t1;

    for (int i = 0; i < 3; i++) {
        if (r1.direction[i] == 0 && (r1.origin[i] < -0.5 || r1.direction[i] > 0.5)) {//parallel to the slab in question
            return pt; //misses the cube
        } else {
            t0 = (-0.5 - r1.origin[i])/ r1.direction[i];
            t1 = (0.5 - r1.origin[i]) / r1.direction[i];
            normNear = glm::vec4(0);
            normNear[i] = -1; //bottom plane (ray coming from bottom)
            normFar = glm::vec4(0);
            normFar[i] = 1;
            if (t0 > t1) {
                //makes t0 represent the intersection w/ slab closer to origin
                float temp = t0;
                t0 = t1;
                t1 = temp;
                normNear[i] = 1;
                normFar[i] = -1;
            }
            if (t0 > t_near) {
                t_near = t0;
                pt.normal = glm::vec3(normNear);
            }
            if (t1 < t_far) {
                t_far = t1;
            }
        }
    }
    if (t_near > t_far) {
        pt.normal = glm::vec3(0); //reset normal
        return pt; // box is missed
    } else {
        // NOW get the actual pt of intersection by plugging in t val:
        if (t_near < 0 && t_far >= 0) { // t_near not valid intersection, use t_far
            glm::vec3 temppt = r1.origin + t_far * r1.direction;
            pt.point = glm::vec3(this->invTransformation * glm::vec4(temppt[0], temppt[1], temppt[2], 1));
            pt.normal = glm::normalize(glm::vec3(this->invTransposeT*normFar));
            pt.t = glm::distance(pt.point,r.origin);
            pt.object_hit = (Geometry*)this;
            return pt;

        } else {
            //make sure the saved pt of intersection is in world view!
            glm::vec3 temppt = r1.origin + t_near * r1.direction;
            pt.point = glm::vec3(this->invTransformation * glm::vec4(temppt[0], temppt[1], temppt[2], 1));
            pt.normal = glm::normalize(glm::vec3(this->invTransposeT*normNear));
            pt.t = glm::distance(pt.point,r.origin);
            pt.object_hit = (Geometry*)this;
            return pt;
        }
    }
}


/**
 * @brief createBoxVertexPositions
 * @brief These are functions that are only defined in this cpp file.
 * @brief They're used for organizational purposeswhen filling the arrays
 * @brief used to hold the vertex and index data.
 */

void BoundingBox::createBoxVertexPositions() {
//    box is sitting on axis: x=out, y=horizontal, z=vertical
    box_vert_pos[0] = glm::vec3(0.5f,-0.5f,-0.5f); //front bottom left
    box_vert_pos[1] = glm::vec3(0.5f,-0.5f,0.5f); //front top left
    box_vert_pos[2] = glm::vec3(0.5f,0.5f,0.5f); //front top right
    box_vert_pos[3] = glm::vec3(0.5f,0.5f,-0.5f); //front bottom right
    box_vert_pos[4] = glm::vec3(-0.5f,-0.5f,-0.5f); //back bottom left
    box_vert_pos[5] = glm::vec3(-0.5f,-0.5f,0.5f); //back top left
    box_vert_pos[6] = glm::vec3(-0.5f,0.5f,0.5f); //back top right
    box_vert_pos[7] = glm::vec3(-0.5f,0.5f,-0.5f); //back bottom right

//    box_vert_pos[0] = glm::vec3(1.0f,-1.0f,-1.0f); //front bottom left
//    box_vert_pos[1] = glm::vec3(1.0f,-1.0f,1.0f); //front top left
//    box_vert_pos[2] = glm::vec3(1.0f,1.0f,1.0f); //front top right
//    box_vert_pos[3] = glm::vec3(1.0f,1.0f,-1.0f); //front bottom right
//    box_vert_pos[4] = glm::vec3(-1.0f,-1.0f,-1.0f); //back bottom left
//    box_vert_pos[5] = glm::vec3(-1.0f,-1.0f,1.0f); //back top left
//    box_vert_pos[6] = glm::vec3(-1.0f,1.0f,1.0f); //back top right
//    box_vert_pos[7] = glm::vec3(-1.0f,1.0f,-1.0f); //back bottom right
}




/**
 * @brief createBoxIndices
 */
void BoundingBox::createBoxIndices () {
    //front face
    box_idx[0] = 0;
    box_idx[1] = 1;

    box_idx[2] = 1;
    box_idx[3] = 2;

    box_idx[4] = 2;
    box_idx[5] = 3;

    box_idx[6] = 3;
    box_idx[7] = 0;

    //back face
    box_idx[8] = 4;
    box_idx[9] = 5;

    box_idx[10] = 5;
    box_idx[11] = 6;

    box_idx[12] = 6;
    box_idx[13] = 7;

    box_idx[14] = 7;
    box_idx[15] = 4;

    //side edges
    box_idx[16] = 1;
    box_idx[17] = 5;

    box_idx[18] = 0;
    box_idx[19] = 4;

    box_idx[20] = 2;
    box_idx[21] = 6;

    box_idx[22] = 3;
    box_idx[23] = 7;
}




void BoundingBox::create() {

    createBoxVertexPositions();
    createBoxIndices();

    for(int i = 0; i < 8; i++){
        box_vert_col[i] = glm::vec3(1.0f,0.0f,0.8f); //fuscia
    }

    count = 24; //index count = 24

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(box_idx, 24 * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(box_vert_pos,8 * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);

    bufCol.allocate(box_vert_col, 8 * sizeof(glm::vec3));
}


GLenum BoundingBox::drawMode() { return GL_LINES;}


bool BoundingBox::isMesh() {
    return false;
}







