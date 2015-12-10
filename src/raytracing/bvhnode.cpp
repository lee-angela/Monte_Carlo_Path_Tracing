#include <raytracing/bvhnode.h>
#include <scene/geometry/BoundingBox.h>
#include <scene/geometry/geometry.h>

BVHnode::BVHnode() {
    this->leftChild = NULL;
    this->rightChild = NULL;
    this->boundingBox = new BoundingBox();
    this->geom = NULL;
}

BVHnode::BVHnode(BVHnode* left, BVHnode* right, BoundingBox* box, Geometry* geom) {
    this->leftChild = left;
    this->rightChild = right;
    this->boundingBox = box;
    this->geom = geom;
}


BVHnode::~BVHnode() {
    if (this->leftChild != NULL) {
        leftChild->~BVHnode();
    }
    if (this->rightChild != NULL) {
        rightChild->~BVHnode();
    }
    this->boundingBox = NULL;
    this->geom = NULL;
}




Intersection BVHnode::getIntersection(Ray r) {
    Intersection pt;
    //Intersection result of bvh node returns Intersection.hit_object = boundingBox
    if (this->boundingBox->getIntersection(r).t >= 0.0f) {
        //if ray hits bounding box for this bvh node, recurse to both children
        Intersection leftPt = this->leftChild->getIntersection(r);
        Intersection rightPt = this->leftChild->getIntersection(r);
        if (leftPt.t > 0 && leftPt.t < rightPt.t) {
            return leftPt;
        } else if (rightPt.t > 0 && rightPt.t < leftPt.t) {
            return rightPt;
        }
        return pt; //if didn't hit left or right child
    }
    return pt;
}


bool BVHnode::hasChildren() {
    if (this->leftChild == NULL && this->rightChild == NULL) {
       return true;
    }
}

bool BVHnode::isLeaf() {
    return (this->leftChild == NULL && this->rightChild == NULL && this->geom != NULL);
}


