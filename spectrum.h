//
// Created by mzy on 2018/11/5.
//

#ifndef RENDERING_SPECTRUM_H
#define RENDERING_SPECTRUM_H

#include "core.h"
#define nSpec 1

class spectrum {
public:
    double c[nSpec];
    spectrum(double v = 0) {
        rep(i, nSpec) c[i] = v;
    }
    spectrum operator + (spectrum const &s) const {
        spectrum ret;
        rep(i, nSpec) ret.c[i] = c[i] + s.c[i];
        return ret;
    }
    spectrum operator * (spectrum const &s) const {
        spectrum ret;
        rep(i, nSpec) ret.c[i] = c[i] * s.c[i];
        return ret;
    }
    spectrum &operator += (spectrum const &s) {
        rep(i, nSpec) c[i] += s.c[i];
        return *this;
    }
    spectrum operator * (double const &d) const {
        spectrum ret;
        rep(i, nSpec) ret.c[i] = c[i] * d;
        return ret;
    }
    spectrum operator / (double const &d) const {
        spectrum ret;
        rep(i, nSpec) ret.c[i] = c[i] / d;
        return ret;
    }
    spectrum &operator /= (double const &d) {
        rep(i, nSpec) c[i] /= d;
        return *this;
    }
    double y() {
        return c[0];
    }
};


#endif //RENDERING_SPECTRUM_H
