//
// Created by mzy on 2018/11/21.
//

#ifndef RENDERING_BSDF_H
#define RENDERING_BSDF_H

#include "spectrum.h"

class BSDF {
public:
    virtual spectrum f(Vec3d wo, Vec3d wi) const = 0;

    virtual spectrum sampleF(Vec3d wo, Vec3d &wi, double &pdf, Vec3d n) const = 0;
};

class lambertianBSDF: public BSDF {
public:
    double R;
    lambertianBSDF(): R(0.3) {}
    spectrum f(Vec3d wo, Vec3d wi) const {
        return spectrum(R / acos(-1));
    }

    spectrum sampleF(Vec3d wo, Vec3d &wi, double &pdf, Vec3d n) const {
        wi = sphereSampler();
        if (wi.dot(n) < 0) wi = -wi;
        pdf = wi.dot(n) / acos(-1);
        return spectrum(R / acos(-1));
    }
};


#endif //RENDERING_BSDF_H
