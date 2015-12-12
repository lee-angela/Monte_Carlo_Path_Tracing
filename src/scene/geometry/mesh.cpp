#include <scene/geometry/mesh.h>
#include <la.h>
#include <tinyobj/tiny_obj_loader.h>
#include <iostream>


float Area(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3)
{
    return glm::length(glm::cross(p1 - p2, p3 - p2)) * 0.5f;
}


void Triangle::ComputeArea()
{
    //Extra credit to implement this
    glm::vec3 point1 = glm::vec3(this->transform.T()*glm::vec4(this->points[0],1.0f));
    glm::vec3 point2 = glm::vec3(this->transform.T()*glm::vec4(this->points[1],1.0f));
    glm::vec3 point3 = glm::vec3(this->transform.T()*glm::vec4(this->points[2],1.0f));

    this->area = Area(point1,point2,point3);
}

void Mesh::ComputeArea()
{
    //Extra credit to implement this
    for (int i = 0; i < this->faces.count(); i++) {
        this->faces[i]->transform = this->transform; //set Triangle's transform to that of the mesh's transform
        this->faces[i]->ComputeArea(); //compute the individual triangle's area
        this->area = this->faces[i]->area + this->area; //combine areas of triangles
    }
}

Intersection Triangle::SamplePoint(float a, float b, Intersection isx) {
    Intersection sp;
    return sp;}


Intersection Mesh::SamplePoint(float a, float b, Intersection isx) {
    Intersection sp;
    return sp;
}


Triangle::Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3):
    Triangle(p1, p2, p3, glm::vec3(1), glm::vec3(1), glm::vec3(1), glm::vec2(0), glm::vec2(0), glm::vec2(0))
{
    for(int i = 0; i < 3; i++)
    {
        normals[i] = plane_normal;
    }
}


Triangle::Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3, const glm::vec3 &n1, const glm::vec3 &n2, const glm::vec3 &n3):
    Triangle(p1, p2, p3, n1, n2, n3, glm::vec2(0), glm::vec2(0), glm::vec2(0))
{}


Triangle::Triangle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3, const glm::vec3 &n1, const glm::vec3 &n2, const glm::vec3 &n3, const glm::vec2 &t1, const glm::vec2 &t2, const glm::vec2 &t3)
{
    plane_normal = glm::normalize(glm::cross(p2 - p1, p3 - p2));
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    normals[0] = n1;
    normals[1] = n2;
    normals[2] = n3;
    uvs[0] = t1;
    uvs[1] = t2;
    uvs[2] = t3;
}


//Returns the interpolation of the triangle's three normals based on the point inside the triangle that is given.
glm::vec3 Triangle::GetNormal(const glm::vec3 &P)
{
    float A = Area(points[0], points[1], points[2]);
    float A0 = Area(points[1], points[2], P);
    float A1 = Area(points[0], points[2], P);
    float A2 = Area(points[0], points[1], P);
    return glm::normalize(normals[0] * A0/A + normals[1] * A1/A + normals[2] * A2/A);
}

glm::vec4 Triangle::GetNormal(const glm::vec4 &position)
{
    glm::vec3 result = GetNormal(glm::vec3(position));
    return glm::vec4(result, 0);
}

glm::vec3 Triangle::ComputeNormal(const glm::vec3 &P)
{}
glm::vec3 Mesh::ComputeNormal(const glm::vec3 &P)
{}

//HAVE THEM IMPLEMENT THIS
//The ray in this function is not transformed because it was *already* transformed in Mesh::GetIntersection
Intersection Triangle::GetIntersection(Ray r){
    //1. Ray-plane intersection
    Intersection result;
    float t =  glm::dot(plane_normal, (points[0] - r.origin)) / glm::dot(plane_normal, r.direction);

    glm::vec3 P = r.origin + t * r.direction;
    //2. Barycentric test
    float S = 0.5f * glm::length(glm::cross(points[1] - points[0], points[2] - points[0]));
    float s1 = 0.5f * glm::length(glm::cross(points[1]-P, points[2]-P))/S;
    float s2 = 0.5f * glm::length(glm::cross(points[2]-P, points[0]-P))/S;
    float s3 = 0.5f * glm::length(glm::cross(points[0]-P, points[1]-P))/S;
    float sum = s1 + s2 + s3;

    if(s1 >= 0 && s1 <= 1 && s2 >= 0 && s2 <= 1 && s3 >= 0 && s3 <= 1){
        result.t = t;
        result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(glm::vec3(P)), material->texture);
        result.object_hit = this;
        result.point = glm::vec3(transform.invT()*glm::vec4(P, 1.0f)); //tranform T to world coords
        result.normal = glm::vec3(glm::normalize(transform.invTransT()*glm::vec4(plane_normal, 1.0f))); //transform normal to world coords
        //TODO: Store the tangent and bitangent
        glm::vec3 changeP1 = this->points[1] - this->points[0];
        glm::vec3 changeP2 = this->points[2] - this->points[0];
        glm::vec2 changeUV1 = this->uvs[1] - this->uvs[0];
        glm::vec2 changeUV2 = this->uvs[2] - this->uvs[0];

        result.tangent = glm::normalize(changeUV2[1]*changeP1 - changeUV1[1]*changeP2);
        result.bitangent =glm::normalize((changeP2 - changeUV2[0]*result.tangent)/changeUV2[1]);
    }
    return result;
}


Intersection Mesh::traverseMeshBVHTree(BVHnode* node, Ray r) {
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
            //tranform all the points
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
            leftpt = traverseMeshBVHTree(node->leftChild,r);
            rightpt = traverseMeshBVHTree(node->rightChild,r);

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
                pt = traverseMeshBVHTree(node->leftChild,r);
            } else if (rightpt.object_hit != NULL){
                pt = traverseMeshBVHTree(node->rightChild,r);
            } else {
                return leftpt;
            }
        }

    } //ray did not hit this boundingbox
    return pt;
}



Intersection Mesh::GetIntersection(Ray r){

    Intersection pt;

    if (this->MeshRootBVHnode != NULL) {
        return traverseMeshBVHTree(this->MeshRootBVHnode, r);
    } else {
        return pt;
    }

//    Ray r_loc = r.GetTransformedCopy(transform.invT());
//    Intersection closest;
//    for(int i = 0; i < faces.size(); i++){
//        Intersection isx = faces[i]->GetIntersection(r_loc);
//        if(isx.object_hit != NULL && isx.t > 0 && (isx.t < closest.t || closest.t < 0)){
//            closest = isx;
//        }
//    }
//    if(closest.object_hit != NULL)
//    {
//        Triangle* tri = (Triangle*)closest.object_hit;
//        glm::vec4 P = glm::vec4(closest.t * r_loc.direction + r_loc.origin, 1);
//        closest.point = glm::vec3(transform.T() * P);
//        closest.normal = glm::normalize(glm::vec3(transform.invTransT() * tri->GetNormal(P)));
//        closest.object_hit = this;
//        closest.t = glm::distance(closest.point, r.origin);//The t used for the closest triangle test was in object space
//        //TODO: Store the tangent and bitangent
//        //tan and bitangent will be already filled by triangle intersections.
//    }
//    return closest;
}

void Mesh::SetMaterial(Material *m)
{
    this->material = m;
    for(Triangle *t : faces)
    {
        t->SetMaterial(m);
    }
}

glm::vec2 Mesh::GetUVCoordinates(const glm::vec3 &point)
{
    return glm::vec2();
}


glm::vec2 Triangle::GetUVCoordinates(const glm::vec3 &point)
{
    float A = Area(points[0], points[1], points[2]);
    float A0 = Area(points[1], points[2], point);
    float A1 = Area(points[0], points[2], point);
    float A2 = Area(points[0], points[1], point);
    return uvs[0] * A0/A + uvs[1] * A1/A + uvs[2] * A2/A;
}

void Mesh::LoadOBJ(const QStringRef &filename, const QStringRef &local_path)
{
    QString filepath = local_path.toString(); filepath.append(filename);
    std::vector<tinyobj::shape_t> shapes; std::vector<tinyobj::material_t> materials;
    std::string errors = tinyobj::LoadObj(shapes, materials, filepath.toStdString().c_str());
    std::cout << errors << std::endl;
    if(errors.size() == 0)
    {
        //Read the information from the vector of shape_ts
        for(unsigned int i = 0; i < shapes.size(); i++)
        {
            std::vector<float> &positions = shapes[i].mesh.positions;
            std::vector<float> &normals = shapes[i].mesh.normals;
            std::vector<float> &uvs = shapes[i].mesh.texcoords;
            std::vector<unsigned int> &indices = shapes[i].mesh.indices;
            for(unsigned int j = 0; j < indices.size(); j += 3)
            {
                glm::vec3 p1(positions[indices[j]*3], positions[indices[j]*3+1], positions[indices[j]*3+2]);
                glm::vec3 p2(positions[indices[j+1]*3], positions[indices[j+1]*3+1], positions[indices[j+1]*3+2]);
                glm::vec3 p3(positions[indices[j+2]*3], positions[indices[j+2]*3+1], positions[indices[j+2]*3+2]);

                Triangle* t = new Triangle(p1, p2, p3);
                if(normals.size() > 0)
                {
                    glm::vec3 n1(normals[indices[j]*3], normals[indices[j]*3+1], normals[indices[j]*3+2]);
                    glm::vec3 n2(normals[indices[j+1]*3], normals[indices[j+1]*3+1], normals[indices[j+1]*3+2]);
                    glm::vec3 n3(normals[indices[j+2]*3], normals[indices[j+2]*3+1], normals[indices[j+2]*3+2]);
                    t->normals[0] = n1;
                    t->normals[1] = n2;
                    t->normals[2] = n3;
                }
                if(uvs.size() > 0)
                {
                    glm::vec2 t1(uvs[indices[j]*2], uvs[indices[j]*2+1]);
                    glm::vec2 t2(uvs[indices[j+1]*2], uvs[indices[j+1]*2+1]);
                    glm::vec2 t3(uvs[indices[j+2]*2], uvs[indices[j+2]*2+1]);
                    t->uvs[0] = t1;
                    t->uvs[1] = t2;
                    t->uvs[2] = t3;
                }
                this->faces.append(t);
            }
        }
        std::cout << "" << std::endl;
        //TODO: .mtl file loading
    }
    else
    {
        //An error loading the OBJ occurred!
        std::cout << errors << std::endl;
    }
}

void Mesh::create(){
    //Count the number of vertices for each face so we can get a count for the entire mesh
        std::vector<glm::vec3> vert_pos;
        std::vector<glm::vec3> vert_nor;
        std::vector<glm::vec3> vert_col;
        std::vector<GLuint> vert_idx;

        unsigned int index = 0;

        for(int i = 0; i < faces.size(); i++){
            Triangle* tri = faces[i];
            vert_pos.push_back(tri->points[0]); vert_nor.push_back(tri->normals[0]); vert_col.push_back(tri->material->base_color);
            vert_pos.push_back(tri->points[1]); vert_nor.push_back(tri->normals[1]); vert_col.push_back(tri->material->base_color);
            vert_pos.push_back(tri->points[2]); vert_nor.push_back(tri->normals[2]); vert_col.push_back(tri->material->base_color);
            vert_idx.push_back(index++);vert_idx.push_back(index++);vert_idx.push_back(index++);
        }

        count = vert_idx.size();
        int vert_count = vert_pos.size();

        bufIdx.create();
        bufIdx.bind();
        bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
        bufIdx.allocate(vert_idx.data(), count * sizeof(GLuint));

        bufPos.create();
        bufPos.bind();
        bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
        bufPos.allocate(vert_pos.data(), vert_count * sizeof(glm::vec3));

        bufCol.create();
        bufCol.bind();
        bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
        bufCol.allocate(vert_col.data(), vert_count * sizeof(glm::vec3));

        bufNor.create();
        bufNor.bind();
        bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
        bufNor.allocate(vert_nor.data(), vert_count * sizeof(glm::vec3));

        MeshRootBVHnode = new BVHnode();

        this->boundingBox = new BoundingBox();
        this->setBoundingBox(); //calls set bounding box for mesh

}

Mesh::~Mesh() {
    this->MeshRootBVHnode->~BVHnode();
}


//This does nothing because individual triangles are not rendered with OpenGL; they are rendered all together in their Mesh.
void Triangle::create(){}





void Triangle::setBoundingBox() {
    //find min/max bounds of this triangle
    glm::vec3 tempPt;

    //this.transform.T() is set by Mesh

    float Xminmax[2];
    float Yminmax[2];
    float Zminmax[2];
    Xminmax[0] = (this->transform.T()*glm::vec4(this->points[0],1))[0]; //set default value as first point
    Xminmax[1] = (this->transform.T()*glm::vec4(this->points[0],1))[0];
    Yminmax[0] = (this->transform.T()*glm::vec4(this->points[0],1))[1];
    Yminmax[1] = (this->transform.T()*glm::vec4(this->points[0],1))[1];
    Zminmax[0] = (this->transform.T()*glm::vec4(this->points[0],1))[2];
    Zminmax[1] = (this->transform.T()*glm::vec4(this->points[0],1))[2];

    for (int i = 0; i < 3; i++) {
        //move all points of triangle according to mesh's tranformation matrix
        tempPt = glm::vec3(this->transform.T()*glm::vec4(this->points[i],1));

        if (tempPt[0] < Xminmax[0]) {
            Xminmax[0] = tempPt[0];
        }
        if (tempPt[0] > Xminmax[1]) {
            Xminmax[1] = tempPt[0];
        }
        if (tempPt[1] < Yminmax[0]) {
            Yminmax[0] = tempPt[1];
        }
        if (tempPt[1] > Yminmax[1]) {
            Yminmax[1] = tempPt[1];
        }
        if (tempPt[2] < Zminmax[0]) {
            Zminmax[0] = tempPt[2];
        }
        if (tempPt[2] > Zminmax[1]) {
            Zminmax[1] = tempPt[2];
        }
    }

    //initialize the boundingBox for this triangle
    this->boundingBox = new BoundingBox();

    this->boundingBox->minBound = glm::vec3(Xminmax[0], Yminmax[0], Zminmax[0]);
    this->boundingBox->maxBound = glm::vec3(Xminmax[1], Yminmax[1], Zminmax[1]);
    this->boundingBox->centerPoint = glm::vec3((this->boundingBox->minBound[0]+this->boundingBox->maxBound[0])/2,
            (this->boundingBox->minBound[1]+this->boundingBox->maxBound[1])/2,
            (this->boundingBox->minBound[2]+this->boundingBox->maxBound[2])/2);

    //find new transformation
    float Xscale = glm::abs(Xminmax[1] - Xminmax[0]); // width of box in any dimension
    float Yscale = glm::abs(Yminmax[1] - Yminmax[0]);
    float Zscale = glm::abs(Zminmax[1] - Zminmax[0]);
    this->boundingBox->setTransformationMats(glm::scale(glm::mat4(1.0f),glm::vec3(Xscale, Yscale, Zscale)), glm::translate(glm::mat4(1.0f), this->boundingBox->centerPoint));
    this->boundingBox->create();
}


void Mesh::setBoundingBox() {
    float Xminmax[2];
    float Yminmax[2];
    float Zminmax[2];
    Xminmax[0] = FLT_MAX;
    Xminmax[1] = FLT_MIN;
    Yminmax[0] = FLT_MAX;
    Yminmax[1] = FLT_MIN;
    Zminmax[0] = FLT_MAX;
    Zminmax[1] = FLT_MIN;
    //set bounding boxes for all the triangle faces within this mesh
    for (int i = 0; i < this->faces.size(); i++) {
        this->faces[i]->transform = this->transform; //set the triangle's tranformation matrices to that of this mesh
        this->faces[i]->setBoundingBox();
        if (this->faces[i]->boundingBox->minBound[0] < Xminmax[0]) {
            Xminmax[0] = this->faces[i]->boundingBox->minBound[0];
        }
        if (this->faces[i]->boundingBox->maxBound[0] > Xminmax[1]) {
            Xminmax[1] = this->faces[i]->boundingBox->maxBound[0];
        }
        if (this->faces[i]->boundingBox->minBound[1] < Yminmax[0]) {
            Yminmax[0] = this->faces[i]->boundingBox->minBound[1];
        }
        if (this->faces[i]->boundingBox->maxBound[1] > Yminmax[1]) {
            Yminmax[1] = this->faces[i]->boundingBox->maxBound[1];
        }
        if (this->faces[i]->boundingBox->minBound[2] < Zminmax[0]) {
            Zminmax[0] = this->faces[i]->boundingBox->minBound[2];
        }
        if (this->faces[i]->boundingBox->maxBound[2] > Zminmax[1]) {
            Zminmax[1] = this->faces[i]->boundingBox->maxBound[2];
        }
    }
    createMeshBVHTree(this->MeshRootBVHnode, this->faces, 0);
    this->boundingBox->minBound = glm::vec3(Xminmax[0], Yminmax[0], Zminmax[0]);
    this->boundingBox->maxBound = glm::vec3(Xminmax[1], Yminmax[1], Zminmax[1]);
    this->boundingBox->centerPoint = glm::vec3((this->boundingBox->minBound[0]+this->boundingBox->maxBound[0])/2,
            (this->boundingBox->minBound[1]+this->boundingBox->maxBound[1])/2,
            (this->boundingBox->minBound[2]+this->boundingBox->maxBound[2])/2);

    //find new transformation
    float Xscale = glm::abs(Xminmax[1] - Xminmax[0]);    // width of box in any dimension
    float Yscale = glm::abs(Yminmax[1] - Yminmax[0]);
    float Zscale = glm::abs(Zminmax[1] - Zminmax[0]);
    this->boundingBox->setTransformationMats(glm::scale(glm::mat4(1.0f),glm::vec3(Xscale, Yscale, Zscale)), glm::translate(glm::mat4(1.0f), this->boundingBox->centerPoint));

    this->boundingBox->create();
}


/**
 * @brief compareXcoords
 * @brief Helper function used for sorting geometry on their bounding box X coords
 * @param a
 * @param b
 * @return -1 if a[x]<=b[x] and 1 if a[x]>b[x]
 */
int compareX(Triangle*& a, Triangle*& b) {
    if (a->boundingBox->maxBound[0] <= b->boundingBox->maxBound[0]) {
        return -1;
    } else {
        return 1;
    }
}

int compareY(Triangle*& a, Triangle*& b) {
    if (a->boundingBox->maxBound[1] <= b->boundingBox->maxBound[1]) {
        return -1;
    } else {
        return 1;
    }
}

int compareZ(Triangle*& a, Triangle*& b) {
    if (a->boundingBox->maxBound[2] <= b->boundingBox->maxBound[2]) {
        return -1;
    } else {
        return 1;
    }
}

//calculates surface areas of all bounding boxes of objs in list
float calcAreas(QList<Triangle*> objs) {
    float totSA = 0.0f;
    for (int i = 0; i < objs.length(); i++) {
        totSA += objs[i]->boundingBox->getSurfArea();
    }
    return totSA;
}

void splitAreas(QList<Triangle*> *firstHalf, QList<Triangle*> *secondHalf, QList<Triangle*> all) {
    int length1 = 1;
    int length2 = all.length()-1;
    int splitIdx = 0; //last index included in the first half of objs

    QList<Triangle*> prev1(all.mid(0,length1));
    QList<Triangle*> prev2(all.mid(1,all.length()-1));

    float prevArea1 = calcAreas(prev1);
    float prevArea2 = calcAreas(prev2);

    for (int i = 0; i < all.length(); i++) {
        length1++; //add one to firstHalf list of objs
        length2--; //remove one from secondHalf list
        QList<Triangle*> temp1(all.mid(0, length1));
        QList<Triangle*> temp2(all.mid(splitIdx, length2));

        float area1 = calcAreas(temp1);
        float area2 = calcAreas(temp2);
        if (area1 > area2) { //when area of first is finally greater than area 2
            //check if the previous split was better balanced:

            float prevDifference = glm::abs(prevArea2 - prevArea1);
            float newDifference = glm::abs(area1 - area2);
            if (prevDifference > newDifference) { //they are NOW more evenly split
                prev1 = temp1;
                prev2 = temp2;
            } else { //they were previously more evenly split
                //current prev1 and prev2 contain ideal split
                break; //dont replace the return value with the new lists
            }
        }
    }
    *firstHalf = prev1;
    *secondHalf = prev2;

}

BVHnode* Mesh::createMeshBVHTree(BVHnode* node, QList<Triangle*> faces, int depth) {
    //create the internal mesh BVH node tree

    //mesh has a separate root node (from scene)
    //if there is only 1 geometry object left in this node, then it is a leaf
    if (faces.count() <= 1) {
        //set the bounding box of the geometry equal to this leaf bvh node
        return new BVHnode(NULL, NULL, faces[0]->boundingBox, faces[0]);

    }

    //get overall boundingBox for all remaining objects in this node
    for (int i = 0; i < faces.count(); i++) {
        //keep adding on bounding boxes of each obj to the overall bbox
        node->boundingBox->combineBoxes(faces[i]->boundingBox);
    }

    //find which axis to divide objects along
    if (glm::mod(depth, 3) == 0) { //multiple of 3 : split along x axis
        //sort according to x coords
        qSort(faces.begin(), faces.end(), compareX);
    } else if (glm::mod(depth, 3) == 1){ //non multiple of 3 : split along y axis
        //sort according to y coords
        qSort(faces.begin(), faces.end(), compareY);
    } else { //glm::mod(depth, 3) == 2 : split along z axis
        //sort according to z coords
        qSort(faces.begin(), faces.end(), compareZ);
    }

    QList<Triangle*> firstHalf;
    QList<Triangle*> secondHalf;
    splitAreas(&firstHalf, &secondHalf, faces);

    node->leftChild = createMeshBVHTree(new BVHnode(), firstHalf, depth+1);
    node->rightChild = createMeshBVHTree(new BVHnode(), secondHalf, depth+1);

    node->boundingBox->create();
    return node;

}

bool Triangle::isMesh() {
    return false;
}

bool Mesh::isMesh() {
    return true;
}

