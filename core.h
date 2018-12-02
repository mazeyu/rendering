//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_CORE_H
#define RENDERING_CORE_H

#include <opencv2/opencv.hpp>
#include <ctime>

using namespace cv;
using namespace std;

#define eps 1e-8
#define oo 1e8
#define rep(i, n) for (int i = 0; i < n; i++)
#define replr(i, l, r) for (int i = l; i <= r; i++)

#define mp make_pair

struct ray {
    Vec3d o, dir;
    ray() {}
    ray(Vec3d o, Vec3d dir): o(o), dir(dir) {}
};

double randomReal();
Vec3d sphereSampler();


#endif //RENDERING_CORE_H