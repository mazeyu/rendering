//
// Created by mzy on 2018/11/21.
//

#ifndef RENDERING_BSDF_H
#define RENDERING_BSDF_H

#include "spectrum.h"

class BSDF {
public:
    spectrum f(Vec3d wo, Vec3d wi) const {
        return spectrum(1);
    }

    spectrum sampleF(Vec3d wo, Vec3d &wi, double &pdf, Vec3d n) {
        wi = sphereSampler();
        pdf = 1.0 / 4 / acos(-1);
        if (wi.dot(n) < 0) wi = -wi;
        return spectrum(0.05);
    }
};
extern BSDF lambert;


#endif //RENDERING_BSDF_H
