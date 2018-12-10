//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_PRIMITIVE_H
#define RENDERING_PRIMITIVE_H

#include "shape.h"



class primitive {

public:
    shape *Shape;
    BSDF *bsdf;

    primitive(shape *Shape): Shape(Shape) {
        bsdf = new lambertianBSDF;
    }
    intersection intersect(ray r) {
        auto ret = Shape->intersect(r);
        ret.bsdf = bsdf;
        return ret;
    }

};


#endif //RENDERING_PRIMITIVE_H
