#ifndef BVHNODE
#define BVHNODE

#include <scene/geometry/BoundingBox.h>
#include <raytracing/intersection.h>
class Geometry;


class BVHnode {
public:
    BVHnode();
    ~BVHnode();
    BVHnode(BVHnode* left, BVHnode* right, BoundingBox* box, Geometry* geom);
    BVHnode* leftChild;
    BVHnode* rightChild;

    BoundingBox* boundingBox; //each bvh node contains instance of boundingbox

    Geometry* geom; //should ONLY be non-NULL when bvh node is a LEAF!

    Intersection getIntersection(Ray r);
    bool hasChildren();
    bool isLeaf();

};

#endif // BVHNODE

