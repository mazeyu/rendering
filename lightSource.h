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


#endif //RENDERING_LIGHTSOURCE_H
