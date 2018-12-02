//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_PRIMITIVE_H
#define RENDERING_PRIMITIVE_H

#include "shape.h"

class primitive {

public:
    shape *Shape;

public:
    primitive(shape *Shape): Shape(Shape) {}

};


#endif //RENDERING_PRIMITIVE_H
