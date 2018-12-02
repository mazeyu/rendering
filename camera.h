//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_CAMERA_H
#define RENDERING_CAMERA_H

#include "core.h"
#include "spectrum.h"

class camera {

protected:
    int h, w;
    double *sampleX, *sampleY;
    double *weight;
    spectrum *s;
    Mat pic;
public:
    ray *rays;
    int nRays;
    spectrum *specs;
    camera() {
        w = 192;
        h = 108;
        sampleX = new double[h * w];
        sampleY = new double[h * w];
        rays = new ray[h * w];
        nRays = h * w;
        specs = new spectrum[h * w];
        s = new spectrum[h * w];
        weight = new double[h * w];
        pic = Mat(h, w, CV_64FC1);
    }
    virtual void genRay() = 0;

    void synthesize(string filename);

};

class pinhole : public camera {
    Vec3d o, axes[3];
public:
    void genRay();

    pinhole(Vec3d o, Vec3d at, Vec3d up):o(o) {


        rep(i, 3) o[i] += (rand() % 10001 - 5000) * eps;

        auto dir = at - o;
        //print(dir);
        dir = dir / norm(dir);
        //print(dir);
        up = up - up.dot(dir) * dir;
        up = up / norm(up);
        auto right = dir.cross(up);
        right = right / norm(right) * w / h;

        //print(dir);
        axes[0] = dir - right - up;
        axes[1] = dir - right + up;
        axes[2] = dir + right - up;
        //print(axes);
    }
};

#endif //RENDERING_CAMERA_H
