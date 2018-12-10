//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_LIGHTSOURCE_H
#define RENDERING_LIGHTSOURCE_H

#include "core.h"
#include "spectrum.h"
#include "shape.h"

class lightSource {
public:

    virtual spectrum Le(ray r) const = 0;

    virtual spectrum sampleLe(ray &r, Vec3d &n, double &pdfDir, double &pdfPos) const = 0;
    virtual spectrum sampleLi(intersection &isect, Vec3d &wi, double &pdf) const = 0;
};

class pointLight: public lightSource {

public:

    Vec3d o;
    spectrum I;
    pointLight(Vec3d o, double I = 1): o(o), I(I) {};

    spectrum Le(ray r) const {
        return spectrum();
    }

    spectrum sampleLe(ray &r, Vec3d &n, double &pdfDir, double &pdfPos) const {
        auto dir = sphereSampler();
        r = ray(o, dir);
        n = dir;
        pdfDir = 1.0 / 4 / acos(-1);
        pdfPos = 1;
        return I;
    }

    spectrum sampleLi(intersection &isect, Vec3d &wi, double &pdf) const {
        wi = o - isect.hit;
        pdf = 1;
        return I / norm(wi) / norm(wi);
    }

};

class areaLight: public lightSource {

public:

    Vec3d o, xAxis, yAxis, n0;
    double R;
    spectrum I;

    areaLight(Vec3d o, Vec3d xAxis, Vec3d yAxis, double R, double I): o(o), xAxis(xAxis), yAxis(yAxis), R(R), I(I) {
        n0 = xAxis.cross(yAxis);
        n0 /= norm(n0);
    };

    spectrum Le(ray r) const {
        return spectrum();
    }

    spectrum sampleLe(ray &r, Vec3d &n, double &pdfDir, double &pdfPos) const {
        auto dir = sphereSampler();
        auto xy = diskSampler();
        r = ray(o + xAxis * xy[0] + yAxis * xy[1], dir);
        n = n0;
        if (dir .dot(n) < 0) dir *= -1;
        pdfDir = dir.dot(n) / acos(-1);
        pdfPos = 1 / (acos(-1) * R * R);
        return I;
    }

    spectrum sampleLi(intersection &isect, Vec3d &wi, double &pdf) const {
        auto xy = diskSampler();
        auto st = o + xAxis * xy[0] + yAxis * xy[1];
        wi = st - isect.hit;
        double normwi = norm(wi);
        pdf = 1 / (cos(-1) * R * R) * normwi * normwi * normwi / abs((isect.nShading.dot(wi))) ;
        if (wi.dot(n0) > 0) return spectrum(0);
        return I;
    }

};


#endif //RENDERING_LIGHTSOURCE_H
