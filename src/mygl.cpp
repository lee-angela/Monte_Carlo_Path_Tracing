#include "mygl.h"
#include <la.h>
#include "raytracing/intersection.h"

#include <iostream>
#include <QApplication>
#include <QKeyEvent>
#include <QXmlStreamReader>
#include <QFileDialog>
#include <renderthread.h>
//#include <tbb/tbb.h>
#include <QElapsedTimer>
#include <raytracing/samplers/stratifiedpixelsampler.h>
#include <scene/geometry/mesh.h>

//using namespace tbb;

MyGL::MyGL(QWidget *parent)
    : GLWidget277(parent)
{
    setFocusPolicy(Qt::ClickFocus);
}

MyGL::~MyGL()
{
    makeCurrent();

    vao.destroy();
}

void MyGL::initializeGL()
{
    // Create an OpenGL context
    initializeOpenGLFunctions();
    // Print out some information about the current OpenGL context
    debugContextVersion();

    // Set a few settings/modes in OpenGL rendering
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    // Set the size with which points should be rendered
    glPointSize(5);
    // Set the color with which the screen is filled at the start of each render call.
    glClearColor(0.5, 0.5, 0.5, 1);

    printGLErrorLog();

    // Create a Vertex Attribute Object
    vao.create();

    // Create and set up the diffuse shader
    prog_lambert.create(":/glsl/lambert.vert.glsl", ":/glsl/lambert.frag.glsl");
    // Create and set up the flat-color shader
    prog_flat.create(":/glsl/flat.vert.glsl", ":/glsl/flat.frag.glsl");

    // We have to have a VAO bound in OpenGL 3.2 Core. But if we're not
    // using multiple VAOs, we can just bind one once.
    vao.bind();

    //Test scene data initialization
    scene.CreateTestScene();
    integrator.scene = &scene;
    integrator.intersection_engine = &intersection_engine;
    intersection_engine.scene = &scene;
    ResizeToSceneCamera();

    //create new tree with this new set of geometry!
    this->intersection_engine.BVHrootNode = new BVHnode();
    this->intersection_engine.BVHrootNode = createBVHtree(this->intersection_engine.BVHrootNode, this->scene.objects, 0);
    update();
}

void MyGL::resizeGL(int w, int h)
{
    gl_camera = Camera(w, h);

    glm::mat4 viewproj = gl_camera.getViewProj();

    // Upload the projection matrix
    prog_lambert.setViewProjMatrix(viewproj);
    prog_flat.setViewProjMatrix(viewproj);

    printGLErrorLog();
}

// This function is called by Qt any time your GL window is supposed to update
// For example, when the function updateGL is called, paintGL is called implicitly.
void MyGL::paintGL()
{
    // Clear the screen so that we only see newly drawn images
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Update the viewproj matrix
    prog_lambert.setViewProjMatrix(gl_camera.getViewProj());
    prog_flat.setViewProjMatrix(gl_camera.getViewProj());
    GLDrawScene();
}

void MyGL::drawTriangles(BVHnode* meshNode, int depth) {

    if ((depth % 5) == 0) {
        prog_flat.setModelMatrix(meshNode->boundingBox->transformation);
        prog_flat.draw(*this, *meshNode->boundingBox);
    }
    //recurse and draw children's boxes
    if(meshNode->leftChild != NULL) {
        drawTriangles(meshNode->leftChild, depth+1);
    }
    if (meshNode->rightChild != NULL) {
        drawTriangles(meshNode->rightChild, depth+1);
    }
}

/**
 * @brief drawBVHnodes
 * @brief Helper function that draws BVHnodes
 */
void MyGL::drawBoxes(BVHnode* node) {
    if (node->geom != NULL && node->geom->isMesh()) {
        drawTriangles(((Mesh*)node->geom)->MeshRootBVHnode, 0);
    }

    prog_flat.setModelMatrix(node->boundingBox->transformation);
    prog_flat.draw(*this, *node->boundingBox);

    //recurse and draw children's boxes
    if(node->leftChild != NULL) {
        drawBoxes(node->leftChild);
    }
    if (node->rightChild != NULL) {
        drawBoxes(node->rightChild);
    }
}

void MyGL::GLDrawScene()
{
    for(Geometry *g : scene.objects)
    {
        if(g->drawMode() == GL_TRIANGLES)
        {
            prog_lambert.setModelMatrix(g->transform.T());
            prog_lambert.draw(*this, *g);
        }
        else if(g->drawMode() == GL_LINES)
        {
            prog_flat.setModelMatrix(g->transform.T());
            prog_flat.draw(*this, *g);
        }
    }
    for(Geometry *l : scene.lights)
    {
        prog_flat.setModelMatrix(l->transform.T());
        prog_flat.draw(*this, *l);
    }
    prog_flat.setModelMatrix(glm::mat4(1.0f));
    prog_flat.draw(*this, scene.camera);

    //Recursively traverse the BVH hierarchy stored in the intersection engine and draw each node
    drawBoxes(this->intersection_engine.BVHrootNode);

}

void MyGL::ResizeToSceneCamera()
{
    this->setFixedWidth(scene.camera.width);
    this->setFixedHeight(scene.camera.height);
    gl_camera = Camera(scene.camera);
}

void MyGL::keyPressEvent(QKeyEvent *e)
{
    float amount = 2.0f;
    if(e->modifiers() & Qt::ShiftModifier){
        amount = 10.0f;
    }
    // http://doc.qt.io/qt-5/qt.html#Key-enum
    if (e->key() == Qt::Key_Escape) {
        QApplication::quit();
    } else if (e->key() == Qt::Key_Right) {
        gl_camera.RotateAboutUp(-amount);
    } else if (e->key() == Qt::Key_Left) {
        gl_camera.RotateAboutUp(amount);
    } else if (e->key() == Qt::Key_Up) {
        gl_camera.RotateAboutRight(-amount);
    } else if (e->key() == Qt::Key_Down) {
        gl_camera.RotateAboutRight(amount);
    } else if (e->key() == Qt::Key_1) {
        gl_camera.fovy += amount;
    } else if (e->key() == Qt::Key_2) {
        gl_camera.fovy -= amount;
    } else if (e->key() == Qt::Key_W) {
        gl_camera.TranslateAlongLook(amount);
    } else if (e->key() == Qt::Key_S) {
        gl_camera.TranslateAlongLook(-amount);
    } else if (e->key() == Qt::Key_D) {
        gl_camera.TranslateAlongRight(amount);
    } else if (e->key() == Qt::Key_A) {
        gl_camera.TranslateAlongRight(-amount);
    } else if (e->key() == Qt::Key_Q) {
        gl_camera.TranslateAlongUp(-amount);
    } else if (e->key() == Qt::Key_E) {
        gl_camera.TranslateAlongUp(amount);
    } else if (e->key() == Qt::Key_F) {
        gl_camera.CopyAttributes(scene.camera);
    } else if (e->key() == Qt::Key_R) {
        scene.camera = Camera(gl_camera);
        scene.camera.recreate();
    }
    gl_camera.RecomputeAttributes();
    update();  // Calls paintGL, among other things
}

void clearTree(BVHnode* root) {
    if (root == NULL) {
        return;
    } else {
        if (root->leftChild != NULL) {
            clearTree(root->leftChild);
        }
        if (root->rightChild != NULL) {
            clearTree(root->rightChild);
        }
    }
    root->~BVHnode();
}


/**
 * @brief compareXcoords
 * @brief Helper function used for sorting geometry on their bounding box X coords
 * @param a
 * @param b
 * @return -1 if a[x]<=b[x] and 1 if a[x]>b[x]
 */
int compareXcoords(const Geometry * a, const Geometry * b) {
    if (a->boundingBox->maxBound[0] <= b->boundingBox->maxBound[0]) {
        return -1;
    } else {
        return 1;
    }
}

int compareYcoords(const Geometry * a, const Geometry * b) {
    if (a->boundingBox->maxBound[1] <= b->boundingBox->maxBound[1]) {
        return -1;
    } else {
        return 1;
    }
}

int compareZcoords(const Geometry* a, const Geometry* b) {
    if (a->boundingBox->maxBound[2] <= b->boundingBox->maxBound[2]) {
        return -1;
    } else {
        return 1;
    }
}


/**
 * @brief MyGL::createBVHtree
 * @brief recursively splits the objects in geometry evenly (along x/y axis) depending on depth
 * @param root
 * @param objs
 * @param depth
 * @return the root BVHnode of the subtree that is being created
 */
BVHnode* MyGL::createBVHtree(BVHnode* root, QList<Geometry*> objs, int depth) {

    //if there is only 1 geometry object left in this node, then it is a leaf
    if (objs.count() <= 1) {
        //set the bounding box of the geometry equal to this leaf bvh node
        return new BVHnode(NULL, NULL, objs[0]->boundingBox, objs[0]);

    }

    //get overall boundingBox for all remaining objects in this node
    for (int i = 0; i < objs.count(); i++) {
        //keep adding on bounding boxes of each obj to the overall bbox
        root->boundingBox->combineBoxes(objs[i]->boundingBox);
    }

    //find which axis to divide objects along
    if (glm::mod(depth, 3) == 0) { //multiple of 3 : split along x axis
        //sort according to x coords
        std::sort(objs.begin(), objs.end(), compareXcoords);
    } else if (glm::mod(depth, 3) == 1){ //non multiple of 3 : split along y axis
        //sort according to y coords
        std::sort(objs.begin(), objs.end(), compareYcoords);
    } else { //glm::mod(depth, 3) == 2 : split along z axis
        //sort according to z coords
        std::sort(objs.begin(), objs.end(), compareZcoords);
    }

    //split the geometry objects according to appropriate coords
    QList<Geometry*> leftHalfObjs(objs.mid(0,floor((float)objs.length()/2.0f)));
    QList<Geometry*> rightHalfObjs(objs.mid(floor((float)objs.length()/2.0f), ceil((float)objs.length()/2.0f)));

    root->leftChild = createBVHtree(new BVHnode(), leftHalfObjs, depth+1);
    root->rightChild = createBVHtree(new BVHnode(), rightHalfObjs, depth+1);

    root->boundingBox->create();
    return root;
}



void MyGL::SceneLoadDialog()
{
    QString filepath = QFileDialog::getOpenFileName(0, QString("Load Scene"), QString("../scene_files"), tr("*.xml"));
    if(filepath.length() == 0)
    {
        return;
    }

    QFile file(filepath);
    int i = filepath.length() - 1;
    while(QString::compare(filepath.at(i), QChar('/')) != 0)
    {
        i--;
    }
    QStringRef local_path = filepath.leftRef(i+1);
    //Reset all of our objects
    scene.Clear();
    integrator = Integrator();

    //delete existing bvh tree
    clearTree(this->intersection_engine.BVHrootNode);
    delete this->intersection_engine.BVHrootNode;

    intersection_engine = IntersectionEngine();

    //Load new objects based on the XML file chosen.
    xml_reader.LoadSceneFromFile(file, local_path, scene, integrator);
    integrator.scene = &scene;
    integrator.intersection_engine = &intersection_engine;
    intersection_engine.scene = &scene;

    //create new tree with this new set of geometry!
    this->intersection_engine.BVHrootNode = new BVHnode();
    this->intersection_engine.BVHrootNode = createBVHtree(this->intersection_engine.BVHrootNode, this->scene.objects, 0);

    update();
}

void MyGL::RaytraceScene()
{
    QString filepath = QFileDialog::getSaveFileName(0, QString("Save Image"), QString("../rendered_images"), tr("*.bmp"));
    if(filepath.length() == 0)
    {
        return;
    }

    glm::vec3 test = integrator.TraceRay(this->scene.camera.Raycast(216, 274), 0);

#define MULTITHREADED
#ifdef MULTITHREADED
    //Set up 16 (max) threads
    unsigned int width = scene.camera.width;
    unsigned int height = scene.camera.height;
    unsigned int x_block_size = (width >= 4 ? width/4 : 1);
    unsigned int y_block_size = (height >= 4 ? height/4 : 1);
    unsigned int x_block_count = width > 4 ? width/x_block_size : 1;
    unsigned int y_block_count = height > 4 ? height/y_block_size : 1;
    if(x_block_count * x_block_size < width) x_block_count++;
    if(y_block_count * y_block_size < height) y_block_count++;

    unsigned int num_render_threads = x_block_count * y_block_count;
    RenderThread **render_threads = new RenderThread*[num_render_threads];

    //start Qtimer

    QElapsedTimer timer;
    timer.start();


    //Launch the render threads we've made
    for(unsigned int Y = 0; Y < y_block_count; Y++)
    {
        //Compute the columns of the image that the thread should render
        unsigned int y_start = Y * y_block_size;
        unsigned int y_end = glm::min((Y + 1) * y_block_size, height);
        for(unsigned int X = 0; X < x_block_count; X++)
        {
            //Compute the rows of the image that the thread should render
            unsigned int x_start = X * x_block_size;
            unsigned int x_end = glm::min((X + 1) * x_block_size, width);
            //Create and run the thread
            render_threads[Y * x_block_count + X] = new RenderThread(x_start, x_end, y_start, y_end, 1, 5, &(scene.film), &(scene.camera), &(integrator));
            render_threads[Y * x_block_count + X]->start();
        }
    }

    bool still_running;
    do
    {
        still_running = false;
        for(unsigned int i = 0; i < num_render_threads; i++)
        {
            if(render_threads[i]->isRunning())
            {
                still_running = true;
                break;
            }
        }
        if(still_running)
        {
            //Free the CPU to let the remaining render threads use it
            QThread::yieldCurrentThread();
        }
    }
    while(still_running);

    //Finally, clean up the render thread objects
    for(unsigned int i = 0; i < num_render_threads; i++)
    {
        delete render_threads[i];
    }
    delete [] render_threads;

#else
    StratifiedPixelSampler pixel_sampler(scene.sqrt_samples);
    for(unsigned int i = 0; i < scene.camera.width; i++)
    {
        for(unsigned int j = 0; j < scene.camera.height; j++)
        {
            QList<glm::vec2> sample_points = pixel_sampler->GetSamples(i, j);
            glm::vec3 accum_color;
            for(int a = 0; a < sample_points.size(); a++)
            {
                glm::vec3 color = integrator.TraceRay(scene.camera.Raycast(sample_points[a]), 0);
                accum_color += color;
            }
            scene.film.pixels[i][j] = accum_color / (float)sample_points.size();
        }
    }
#endif

    //stop Qtimer
    qDebug()<<"this raytrace took "<<timer.elapsed()<<"milliseconds";


    scene.film.WriteImage(filepath);
}
