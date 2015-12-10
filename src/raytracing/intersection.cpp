#include <raytracing/intersection.h>
#include <raytracing/intersectionengine.h>

Intersection::Intersection():
    point(glm::vec3(0)),
    normal(glm::vec3(0)),
    tangent(glm::vec3(0)),
    bitangent(glm::vec3(0)),
    t(-1),
    texture_color(glm::vec3(1.0f))
{
    object_hit = NULL;
}

IntersectionEngine::IntersectionEngine()
{
    scene = NULL;

    //initialize the root node of bvh tree
    BVHrootNode = new BVHnode();
}


/**
 * Helper function that recursively traverses the BVH tree
 *
 * input: BVHnode* node, Ray r
 * output: Intersection pt,
 *          (NOTE: this pt can be default if input node is NULL or there is no found intersection,
 *                  object_hit == NULL, and t = -1 for default pt)
 */

Intersection traverseBVH(BVHnode* node, Ray r) {
    Intersection pt;
    Intersection nothing;

    if (node == NULL) {
        return pt;
    }
    //check for intersection with this node's boundingBox
    pt = node->boundingBox->getIntersection(r);

    // if ray hits this bounding box
    if (pt.object_hit != NULL) { // ELSE pt.object_hit = NULL.
        //check if this is leaf node
        if (node->isLeaf()) {
            return node->geom->GetIntersection(r);
        }

        Intersection leftpt = node->leftChild->boundingBox->getIntersection(r);
        Intersection rightpt = node->rightChild->boundingBox->getIntersection(r);
        //if not leaf node, recurse down to children

        //hits neither of children
        if (leftpt.object_hit == NULL && rightpt.object_hit == NULL) {
            return nothing;
        }
        //IF BOTH CHILDREN ARE INTERSECTED BY RAY, need to check both
        if (leftpt.object_hit != NULL && rightpt.object_hit != NULL) {
            leftpt = traverseBVH(node->leftChild,r);
            rightpt = traverseBVH(node->rightChild,r);

            if (leftpt.t > 0 && rightpt.t < 0) { pt = leftpt; } //left child's geom intersection is closer

            else if (rightpt.t > 0 && leftpt.t < 0) { pt = rightpt;} //right child's geom intersection is closer

            else if (leftpt.t > 0 && rightpt.t > 0) { //ray hit both children's geom
                //take closer geom intersection
                if (leftpt.t < rightpt.t) { pt = leftpt;}
                else { pt = rightpt;}
            } else {
                return nothing;
            }
        }
        //JUST ONE (or none) CHILD IS INTERSECTED
        else {
            if (leftpt.object_hit != NULL) {
                pt = traverseBVH(node->leftChild,r);
            } else if (rightpt.object_hit != NULL){
                pt = traverseBVH(node->rightChild,r);
            } else {
                return leftpt;
            }
        }

    } //ray did not hit this boundingbox
    return pt;
}



Intersection IntersectionEngine::GetIntersection(Ray r)
{

    Intersection pt;

    if (this->BVHrootNode != NULL) {
        return traverseBVH(BVHrootNode, r);
    } else {
        return pt;
    }

// implementation of
//    Intersection nearest;
//    for(Geometry* g : scene->objects)
//    {
//        Intersection isx = g->GetIntersection(r);
//        if((isx.t < nearest.t && isx.object_hit != NULL) || nearest.t < 0)
//        {
//            nearest = isx;
//        }
//    }
//    return nearest;
}

bool IntersectionComp(const Intersection &lhs, const Intersection &rhs)
{
    return lhs.t < rhs.t;
}

QList<Intersection> IntersectionEngine::GetAllIntersections(Ray r)
{
    QList<Intersection> result;
    for(Geometry* g : scene->objects)
    {
        Intersection isx = g->GetIntersection(r);
        if(isx.t > 0)
        {
            result.append(isx);
        }
    }
    std::sort(result.begin(), result.end(), IntersectionComp);
    return result;
}
