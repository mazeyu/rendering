//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_SCENE_H
#define RENDERING_SCENE_H
#include "core.h"
#include "primitive.h"
#include "lightSource.h"
#include "camera.h"
#include "spectrum.h"
#include "thread"

struct bb {
    Vec3d min, max;
    bb(Vec3d min, Vec3d max): min(min), max(max) {}
    //bb(const bb& BB): min(BB.min), max(BB.max) {}
};

struct kdNode {
    kdNode *lc, *rc;
    bb BB;
    vector<primitive*> prims;
    kdNode(bb BB, vector<primitive*> prims): BB(BB), prims(prims) {}
};

class scene {
    int kdCnt = 0;
    vector<primitive*> objs;
    vector<lightSource*> lights;
    camera *cam;
    spectrum integL(ray r);
    spectrum RT(ray r);
    kdNode *kdRoot;
    intersection intersect(ray r);
    vector<primitive*> kdIntersect(kdNode *cur, ray r);
    spectrum UniformSampleOneLight(intersection &isect);

public:
    void init();
    void render();
    void SPPM();
    void addFromFile();
    void addBezier();
    void buildKdTree();

};


#endif //RENDERING_SCENE_H
